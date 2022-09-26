using StatsBase
using ImageMorphology
using LsqFit

WaveletMatrix = dimmatrix(Ti, :𝑓) # Type for DimArrays containing wavelet transform info
export WaveletMatrix

abstract type AbstractBurst end
BurstVector = AbstractVector{<:AbstractBurst}

Base.@kwdef mutable struct Burst <: AbstractBurst
    mask::WaveletMatrix
    thresh
    peak = nothing
    fit = nothing
    width = nothing
    significance = nothing
end

Burst(mask, thresh; kwargs...) = Burst(; mask, thresh, kwargs...)
Burst(mask, thresh, peak; kwargs...) = Burst(; mask, thresh, peak, kwargs...)

duration(B::Burst) = B.fit.param[4]
spectralwidth(B::Burst) = B.fit.param[5]
binarymask(B::Burst) = B.mask .> B.thresh[1]

function threshold(res, thresh, mode)
    if mode == :percentile
        cutoff = percentile(res[:], thresh)
    elseif mode == :std
        cutoff = mean(res[:]) + thresh*std(res[:])
    end
end

"""
Threshold a wavelet spectrum using either a percentile cutoff (`mode=:percentile`) or a standard deviation cutoff (`mode=:std`) of either each frequency band (`eachfreq=true`) or the entire spectrum.
"""
function burstthreshold!(res::WaveletMatrix, thresh; mode=:percentile, eachfreq=true)
    if eachfreq
        cutoffs = threshold.(eachcol(res), thresh, mode)
    else
        cutoffs = threshold(res, thresh, mode)
    end
    res[res .< cutoffs'] .= 0.0
end

burstthreshold(res::WaveletMatrix, thresh; kwargs...) = (y = deepcopy(res); burstthreshold!(y, thresh; kwargs...); y)

function _detectbursts(_res::WaveletMatrix; areacutoff=1)
# Get the connected component masks
    components = ImageMorphology.label_components(_res|>Array)
    areas = ImageMorphology.component_lengths(components)
    idxs = areas .≥ areacutoff
    centroids = ImageMorphology.component_centroids(components)
    centroids = [round.(Int, c) for c in centroids]
    𝑓 = dims(res, Dim{:𝑓})
    t = dims(res, Ti)
    peaks = [(t[c[1]], 𝑓[c[2]]) for c in centroids]
    peaks = [(p[1], p[2], res[Ti(At(p[1])), Dim{:𝑓}(At(p[2]))]) for p in peaks]
    masks = burstmask.((res,), peaks)
    return Burst.(masks, ((minimum(res[_res]), mode, thresh),), peaks)[idxs]
end

"""
Detect bursts from a supplied wavelet spectrum, using thresholding
"""
function detectbursts(res::WaveletMatrix, thresh=99.5; mode=:percentile, kwargs...)
    _res = burstthreshold(res, thresh; mode) .> 0
    _detectbursts(res; kwargs...)
end

"""
This is the preferred, surrogate-based method of detecting bursts. Assumes the surrogate will be stationary, so averages statistics over a wavelet transform of a single surrogate.
`thresh` is in SDs
"""
function detectbursts(x::LFPVector; pass=[30, 100], areacutoff=1, thresh=2, kwargs...) # surrodur=min(length(x), round(Int, 50/step(dims(x, Ti))/minimum(pass))), N=100
    s = surrogate(x, IAAFT())
    # S = [s() for _ in 1:N]
    γₛ = gammafilter(s; pass)
    res = wavelettransform(γₛ)
    # agg = [getindex.(res, (i,)) for i in CartesianIndices(res[1])]
    # agg = DimArray(agg, dims(res[1]))
    σ = mapslices(std, res, dims=1)
    μ = mapslices(mean, res, dims=1)

    # * Determine a threshold from the surrogate wavelet spectra. Lets just use SD's for simplicity
    # * This can be done by normalizing each frequency band with the mean and SD from the surrogate
    # .................will need to At() the wavelet values.............

    γ = gammafilter(x; pass)
    res = wavelettransform(γ; kwargs...)
    _detectbursts(res; areacutoff)
end

function fit!(B::Burst)
    mask = B.mask
    B.fit = fitdiagonalgaussian(mask)
end

function fitdiagonalgaussian(x, y, Z::AbstractMatrix)
    x, y, Z = collect.((x, y, Z))
    xy = hcat(collect.(Iterators.product(x, y))...)'
    z = Z[:]
    function gaussian2(xy, p)
        x, y = eachcol(xy)
        A, μ₁, μ₂, σ₁, σ₂ = p
        y = A.*exp.(-0.5.*(((x .- μ₁)./σ₁).^2 .+ ((y .- μ₂)./σ₂).^2))
    end
    p0 = [  maximum(Z), # A
            mean(x), # μ₁
            mean(y), # μ₂
            (x |> diff |> collect |> extrema |> first)/2, # σ₁
            (y |> diff |> collect |> extrema |> first)/2] # σ₂
    fit = LsqFit.curve_fit(gaussian2, xy, z, p0)
end

function fitdiagonalgaussian(mask::AbstractDimArray)
    𝑓 = dims(mask, 2)
    t = dims(mask, 1)
    f = fitdiagonalgaussian(t, 𝑓, mask|>Array)
end
















# """
# Get a sufficient image of a burst from an entire wavelet spectrum. A threshold at which the power is 1% of the peak power should suffice, but really the threshold should just be small enough for the exact value not to matter.
# Inputs are `res`, the complete wavelet spectrum, `peak`, the (time, freq, value) or (time, freq) tuple of the peak, and `thresh`, the proportion threshold on the peak power. `diffthresh` is a threshold on the gradient of the spectrum, to avoid capturing a bimodal burst.

# Could consider doing this in two passess, to be computationally efficient: first pass by looking along perpendicular directions, with a low threshold, then second pass by fitting a Gaussian on this loose mask and setting a higher threshold.
# """
function burstmask(res::WaveletMatrix, peak; thresh=0.8, n=3)#, diffthresh=0.01)
    pidx = ((Ti∘At)(peak[1]), (Dim{:𝑓}∘At)(peak[2]))
    A = res[pidx...]
    thresh = thresh*A
    # Find the rectangular boundary that contains the threshold. Work in each dimension from the centre. Assumes the burst profile is blobby around the peak and its horizontal/vertical axes.
    𝑓 = dims(res, Dim{:𝑓})
    t = dims(res, Ti)

    bounds = zeros(4) # (tmin, tmax, fmin, fmax)

    # ! This method is not that great
    # ups = res[pidx[2]][t .< peak[1]]
    # dups = -vcat(diff(ups), [0])
    # up = findlast((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[1] = isnothing(up) ? 1 : up

    # ups = res[pidx[2]][t .> peak[1]]
    # dups = vcat([0], diff(ups))
    # up = findfirst((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[2] = isnothing(up) ? size(res, Ti) : findfirst(t .> peak[1]) + up - 1

    # ups = res[pidx[1]][𝑓 .< peak[2]]
    # dups = -vcat(diff(ups), [0])
    # up = findlast((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[3] = isnothing(up) ? 1 : up

    # ups = res[pidx[1]][𝑓 .> peak[2]]
    # dups = vcat([0], diff(ups))
    # up = findfirst((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[4] = isnothing(up) ? size(res, Dim{:𝑓}) : findfirst(𝑓 .> peak[2]) + up - 1

    # bounds[1:2] .= t[Int.(bounds[1:2])]
    # bounds[3:4] .= 𝑓[Int.(bounds[3:4])]

    # * This one looks for the hwhm in each direction, then just multiplies that by a few for good coverage. It's not a terrible thing if we are too broad.
    p = t .< peak[1]
    ups = res[pidx[2]][p]
    up = findlast((ups .< thresh))
    w = length(ups) - up
    up = findlast(p) - n*w
    bounds[1] = max(1, up)

    p = t .> peak[1]
    ups = res[pidx[2]][p]
    up = findfirst((ups .< thresh))
    w = up
    up = findfirst(p) + n*w
    bounds[2] = min(size(res, Ti), up)

    p = 𝑓 .< peak[2]
    ups = res[pidx[1]][p]
    up = findlast((ups .< thresh))
    w = length(ups) - up
    up = findlast(p) - n*w
    bounds[3] = max(1, up)

    p = 𝑓 .> peak[2]
    ups = res[pidx[1]][p]
    up = findfirst((ups .< thresh))
    w =  up
    up = findfirst(p) + n*w
    bounds[4] = min(size(res, Dim{:𝑓}), up)

    bounds[1:2] .= t[Int.(bounds[1:2])]
    bounds[3:4] .= 𝑓[Int.(bounds[3:4])]

    mask = res[Ti(Interval(bounds[1:2]...)), Dim{:𝑓}(Interval(bounds[3:4]...))]
end
