using StatsBase
using ImageMorphology
using LsqFit
using LinearAlgebra
using ProgressLogging
using Interpolations
using HypothesisTests
# using FiniteDifferences

abstract type AbstractBurst end
BurstVector = AbstractVector{<:AbstractBurst}

Base.@kwdef mutable struct Burst <: AbstractBurst
    mask::LogWaveletMatrix
    thresh
    peak = nothing
    fit = nothing
    width = nothing
    significance = nothing
end

Burst(mask, thresh; kwargs...) = Burst(; mask, thresh, kwargs...)
Burst(mask, thresh, peak; kwargs...) = Burst(; mask, thresh, peak, kwargs...)

duration(B::Burst) = B.fit.param[4]
logspectralwidth(B::Burst) = B.fit.param[5]
function spectralwidth(B::Burst)
    σ = logspectralwidth(B)
    μ = logpeakfreq(B)
    return (exp10(σ^2) - 1)*exp10(2*μ + σ^2)
end
amplitude(B::Burst) = B.fit.param[1]
mask(B::Burst) = B.mask
# binarymask(B::Burst) = B.mask .> B.thresh[1]
peaktime(B::Burst) = B.fit.param[2]
maxfreq(B::Burst) = exp10(dims(mask(B), 2)[findmax(mask(B))[2][2]])
maxtime(B::Burst) = dims(mask(B), 1)[findmax(mask(B))[2][1]]
logpeakfreq(B::Burst) = B.fit.param[3]
peakfreq(B::Burst) = B |> logpeakfreq |> exp10
df(B::Burst) = mean(diff(collect(dims(mask(B), Dim{:logfrequency}))))
dt(B::Burst) = step(dims(mask(B), Ti))
maskduration(B::Burst) = dims(mask(B), Ti) |> extrema |> collect |> diff |> first
maskspectralwidth(B::Burst) = dims(mask(B), Dim{:logfrequency}) |> extrema |> collect |> diff |> first
fiterror(B::Burst) = std(B.fit.resid./B.mask[:])
interval(B::Burst) = peaktime(B)±(0.5*duration(B))
inany(x, V::Vector{<:AbstractInterval}) = any(in.((x,), V))

function basicfilter!(B::BurstVector; fmin=1, tmin=0.008, tmax=1, pass=[0, Inf])
    fmin = log10(fmin)
    # pass = Interval(pass...)
    # passes(x) = all([y in pass for y in x])
    filter!(b->size(mask(b), Ti) > tmin/dt(b), B)
    filter!(b->size(mask(b), Ti) < tmax/dt(b), B)
    filter!(b->size(mask(b), Dim{:logfrequency}) > fmin/df(b), B)
    # filter!(b->passes(extrema(dims(mask(b), Dim{:frequency}))), B)
end

function bandfilter!(B::BurstVector; pass=[50, 60])
    pass = Interval(pass...)
    filter!(b->maxfreq(b) ∈ pass, B)
end

bandfilter!(; kwargs...) = x->bandfilter!(x; kwargs...)

function filtergammabursts!(B::BurstVector; fmin=1, tmin=0.008, pass=[50, 60]) # fmin in Hz, tmin in s
    fmin = log10(fmin)
    pass = Interval(pass...)
    filter!(b->size(mask(b), Ti) > tmin/dt(b), B)
    filter!(b->size(mask(b), Dim{:logfrequency}) > fmin/df(b), B)
    filter!(b->peakfreq(b) ∈ pass, B)
    # ..... Other stuff
end

# # Old method, global distribution
# function threshold(res, thresh, method)
#     if method == :percentile
#         cutoff = percentile(res[:], thresh)
#     elseif method == :std
#         cutoff = thresh*std(res[:])
#     end
# end
function threshold(res, thresh, method)
    if method == :percentile
        cutoff = mapslices(x->percentile(x, thresh), res, dims=Ti)
    elseif method == :std
        cutoff = thresh*mapslices(std, res, dims=Ti)
    elseif method == :iqr
        cutoff = thresh*mapslices(iqr, res, dims=Ti)./1.35 # 1.35 to be consistent with std threshold on normal data
    end
end

function surrogatethreshold()
end

# grad(x::Real, y::Real) = x - y
# grad(x::AbstractVector, h) = grad.((@view x[3:end]), (@view x[1:end-2]))./(2*h)
# grad(h::Real) = x -> grad(x, h)

"""
Threshold a wavelet spectrum using either a percentile cutoff (`method=:percentile`) or a standard deviation cutoff (`method=:std`) of either each frequency band (`eachfreq=true`) or the entire spectrum. You probably want to FOOOF the spectrum before this
"""
function burstthreshold!(res::LogWaveletMatrix, thresh; method=:iqr, zerograd=0.0)
    @assert dims(res, Ti).val.data isa AbstractRange "Rectify the LFP array before calculating the wavelet transform"
    cutoffs = threshold(res, thresh, method)
    res[res .< cutoffs] .= 0.0
end

burstthreshold(res, thresh; kwargs...) = (y = deepcopy(res); burstthreshold!(y, thresh; kwargs...); y)

burstthreshold!(res::WaveletMatrix, thresh; kwargs...) = burstthreshold!(convert(LogWaveletMatrix, res), thresh; kwargs...)


"""

`thresh` is a proportion of the average gradient of the wavelet spectrum below which a gradie
"""
function burstcurvature!(res::LogWaveletMatrix, thresh=0)
    @assert dims(res, Ti).val.data isa AbstractRange "Rectify the LFP array before calculating the wavelet transform"

    xs = dims(res, Ti)
    @assert xs.val.data isa AbstractRange "Expected rectified time indices"
    xs = xs.val.data
    ys = collect(dims(res, Dim{:logfrequency}))
    @assert std(diff(ys))/std(ys) < 1e-3 "Logarithmic frequency bins are not approximately uniform"
    df = median(diff(ys))
    _ys = minimum(ys):df:(maximum(ys)+3.0*df)
    ys = _ys[1:length(ys)]
    A = collect(res)
    nodes = (xs, ys)
    itp = scale(interpolate(A, BSpline(Quadratic(Line(OnGrid())))), nodes...)
    grad = (x, y) -> Interpolations.gradient(itp, x, y)
    hess = (x, y) -> Interpolations.hessian(itp, x, y)
    # G = [grad(x, y) for x in xs, y in ys]
    H = [hess(x, y) for x in xs, y in ys]
    kx = [h[1, 1] for h in H] # Curvature in time
    ky = [h[2, 2] for h in H] # Curvature in frequency
    # N = norm.(G)
    # N̄ = mean(N)
    res[kx .> -thresh.*median(kx)] .= 0.0
    res[ky .> -thresh.*median(ky)] .= 0.0
end
burstcurvature(res, thresh=0; kwargs...) = (y = deepcopy(res); burstcurvature!(y, thresh; kwargs...); y)
burstcurvature!(res::WaveletMatrix, thresh=0; kwargs...) = burstcurvature!(convert(LogWaveletMatrix, res), thresh; kwargs...)



function widen(x, δ=0.5; upperbound=[Inf, Inf])
    δ = δ/2
    @assert length(x[1]) == length(x[2]) == 2
    Δ = [x[2][1] - x[1][1], x[2][2] - x[1][2]]
    return [max.(1, floor.(Int, x[1] .- δ.*Δ)), min.(upperbound, ceil.(Int, x[2] .+ δ.*Δ))]
end

function _detectbursts(res::LogWaveletMatrix; thresh=3, curvaturethresh=1, boundingstretch=0.5, method=:iqr, areacutoff=1)
    # @info "Thresholding amplitudes"
    _res = burstthreshold(res, thresh; method) .> 0
    # @info "Thresholding curvatures"
    _res = _res .& (burstcurvature(res, curvaturethresh) .> 0)

    # @info "Finding connected components"
    components = ImageMorphology.label_components(_res|>Array)
    areas = ImageMorphology.component_lengths(components)[2:end]
    idxs = areas .≥ areacutoff
    centroids = ImageMorphology.component_centroids(components)[2:end]
    centroids = [round.(Int, c) for c in centroids]
    𝑓 = dims(_res, Dim{:logfrequency})
    t = dims(_res, Ti)
    peaks = [(t[c[1]], 𝑓[c[2]]) for c in centroids]
    peaks = [(p[1], p[2], _res[Ti(At(p[1])), Dim{:logfrequency}(At(p[2]))]) for p in peaks]
    # ! Might want to look for gaussian peak instead
    # masks = burstmask.((_res,), peaks)
    bb = ImageMorphology.component_boxes(components)[2:end]
    bb = [widen(b, boundingstretch; upperbound=size(res)) for b in bb]
    masks = [res[b[1][1]:b[2][1], b[1][2]:b[2][2]] for b in bb]
    B = Burst.(masks, ((method, thresh, curvaturethresh),), peaks)[idxs]
end

function mmap_detectbursts(res::LogWaveletMatrix; window=50000, kwargs...)
    ti = _slidingwindow(collect(dims(res, Ti)), window; tail=true)
    ti = [Interval(extrema(t)...) for t in ti]
    B = Vector{Burst}()
    threadlog, threadmax = (0, length(ti))
    @withprogress name="Burst detection" begin
        for ts in ti
            subres = res[Ti(ts)]
            append!(B, _detectbursts(subres; kwargs...))
            if threadmax > 1
                Threads.threadid() == 1 && (threadlog += 1)%1 == 0 && @logprogress threadlog/threadmax
            end
        end
    end
    return B
end

function pvalues(B::BurstVector, Bs::BurstVector, f::Function; test=OneSampleTTest, tail=:left)
    d = f.(B)
    ds = f.(Bs)
    p = pvalue.(test.((ds,), d); tail).*length(d) # Bonferroni correction
end

function pvalues(B::BurstVector, Bs::BurstVector, fs::AbstractVector{<:Function}; test=OneSampleHotellingT2Test, tail=nothing)
    d = collect.(zip([f.(B) for f in fs]...))
    ds = hcat([f.(Bs) for f in fs]...)
    tests = test.((ds,), d)
    p = pvalue.(tests).*length(d) # Bonferroni correction
end

function significancefilter!(B::BurstVector, Bs::BurstVector, f; thresh=0.05, kwargs...)
    p = pvalues(B, Bs, f; kwargs...)
    deleteat!(B, p .> thresh)
end

"""
Filter a vector of bursts based of a vector of surrogate bursts
"""
function surrogatefilter!(B::BurstVector, Bs::BurstVector; stats = [duration, peakfreq], joint=true)
    if joint
        significancefilter!(B, Bs, stats; test=OneSampleHotellingT2Test, tail=nothing)
    else
        for s in stats
            # ! The problem is, the distribution changes with each of these. Will want to take the maximum p-value, rather successively filtering
            significancefilter!(B, Bs, s; test=OneSampleTTest, tail=:left, thresh=0.01)
        end
    end
end

"""
Detect bursts from a supplied wavelet spectrum, using thresholding
`boundingstretch` increases the bounding box slightly so for a more accurate fit. Give as a proportion of the threshold bounding box
`detection` can be `_detectbursts` or `mmap_detectbursts`
"""
function detectbursts(res::LogWaveletMatrix; pass=nothing, dofit=true, detection=_detectbursts, kwargs...)
    B = detection(res; kwargs...)
    basicfilter!(B)
    isnothing(pass) || (@info "Filtering in the $(pass) Hz band"; bandfilter!(B; pass))

    if dofit
        @info "Fitting burst profiles"
        fit!(B)
        sort!(B, by=peaktime)
    end
    return B
end

function detectbursts(x::LFPVector; kwargs...) # surrodur=min(length(x), round(Int, 50/step(dims(x, Ti))/minimum(pass))), N=100
    # s = surrogate(x, AP())
    # # S = [s() for _ in 1:N]
    # # γₛ = gammafilter(s; pass)
    # res = wavelettransform(γₛ)
    # agg = [getindex.(res, (i,)) for i in CartesianIndices(res[1])]
    # agg = DimArray(agg, dims(res[1]))
    # σ = mapslices(std, res, dims=1)
    # μ = mapslices(mean, res, dims=1)

    # γ = gammafilter(x; pass)
    res = fooofedwavelet(x)
    detectbursts(res; kwargs...)
end

function fit!(B::Burst)
    mask = B.mask

    # The fitting takes a while. Let's downsample a little bit if we can
    while size(mask, 1) > 100
        mask = mask[1:2:end, :]
    end
    while size(mask, 2) > 100
        mask = mask[:, 1:2:end]
    end
    B.fit = fitdiagonalgaussian(mask)
end

function fit!(ℬ::AbstractVector{<:AbstractBurst})
    threadlog, threadmax = (0, length(ℬ)/Threads.nthreads())
    @withprogress name="Fitting bursts" begin
        Threads.@threads for B in ℬ
            try
                fit!(B)
            catch e
                @warn "Could not fit a burst. Skipping"
            end
            Threads.threadid() == 1 && (threadlog += 1)%1 == 0 && @logprogress threadlog/threadmax
        end
    end
    deleteat!(ℬ, findall(.!isa.(getfield.(ℬ, (:fit,)), (LsqFit.LsqFitResult,))))
end

function gaussian2!(F, xy, p)
    x, y = eachcol(xy)
    A, μ₁, μ₂, σ₁, σ₂ = p
    @. F = A*exp(-0.5*(((x - μ₁)/σ₁)^2 + ((y - μ₂)/σ₂)^2))
end
gaussian2(xy::Tuple, p) = (F = [0.0]; gaussian2!(F, collect(xy)', p) |> first)

function gaussian2_j!(J::Array{Float64,2}, xy, p)
    x, y = eachcol(xy)
    A, μ₁, μ₂, σ₁, σ₂ = p
    @. J[:,1] .= exp(-0.5*(((x - μ₁)/σ₁)^2 + ((y - μ₂)/σ₂)^2))                       #dF/A
    @. @views J[:,2] .= A*(x-μ₁)*J[:,1]/(σ₁^2)         #dF/μ₁
    @. @views J[:,3] .= A*(y-μ₂)*J[:,1]/(σ₂^2)         #dF/μ₂
    @. @views J[:,4] .= A*(x-μ₁)^2*J[:,1]/(σ₁^3)       #dF/σ₁
    @. @views J[:,5] .= A*(x-μ₂)^2*J[:,1]/(σ₂^3)       #dF/σ₂
end

function fitdiagonalgaussian(x, y, Z::AbstractMatrix)
    x, y = collect.((x, y))
    xy = hcat(collect.(Iterators.product(x, y))...)'
    z = Z[:]
    A0, μ0 = findmax(z)
    p0 = [  A0, # A
            xy[μ0, 1], # μ₁
            xy[μ0, 2], # μ₂
            (x |> extrema |> collect |> diff |> first)/3, # σ₁
            (y |> extrema |> collect |> diff |> first)/3] # σ₂
    T = eltype(z)
    xy = xy .|> T
    p0 = p0 .|> T
    fit = LsqFit.curve_fit(gaussian2!, xy, z, p0; inplace=true, autodiff=:forwarddiff)
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
function burstmask(res::LogWaveletMatrix, peak; thresh=0.8, n=3)#, diffthresh=0.01)
    pidx = ((Ti∘At)(peak[1]), (Dim{:logfrequency}∘At)(peak[2]))
    A = res[pidx...]
    thresh = thresh*A
    # Find the rectangular boundary that contains the threshold. Work in each dimension from the centre. Assumes the burst profile is blobby around the peak and its horizontal/vertical axes.
    𝑓 = dims(res, Dim{:logfrequency})
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
    # bounds[4] = isnothing(up) ? size(res, Dim{:frequency}) : findfirst(𝑓 .> peak[2]) + up - 1

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
    bounds[4] = min(size(res, Dim{:frequency}), up)

    bounds[1:2] .= t[Int.(bounds[1:2])]
    bounds[3:4] .= 𝑓[Int.(bounds[3:4])]

    mask = res[Ti(Interval(bounds[1:2]...)), Dim{:frequency}(Interval(bounds[3:4]...))]
end


"""
Pair each burst in A with a burst in B, i.e. where the burst in A is closest in time to the burst in B. Not all bursts will be paired.
"""
function pairbursts(A::BurstVector, B::BurstVector; feed=:forward)
    B = sort(deepcopy(B), by=peaktime)
    A = sort(deepcopy(A), by=peaktime)
    C = []
    D = []
    ta = peaktime.(A)
    tb = peaktime.(B)
    for (i, a) in enumerate(A)
        ts = ta[i] .> tb
        t = findlast(ts)
        if ~isnothing(t)
            if feed == :forward
                ts[t+1] = true
                t = findlast(ts)
            end
            push!(C, a)
            push!(D, B[t])
            deleteat!(B, ts) # Remove all bursts before `a``
            deleteat!(tb, ts)
        end
    end
    return C, D
end


function burstdelta(A::BurstVector, B::BurstVector; feed=:forward)
    C, D = pairbursts(A, B; feed)
    Δ = peaktime.(D) .- peaktime.(C)
end


function burstspikestats(B, Sp, channels; sessionid, probeid, kwargs...)
    unitchannels = AN.getclosestchannels(sessionid, probeid, keys(Sp), channels)
    # * Count the number of spikes within bursts vs outside bursts. Assume Sp contains spikes only in the duration of the LFP used to calculate bursts.
    cols = [:unit, :channel, :mean_rate, :burst_rate, :nonburst_rate]
    stats = DataFrame([[] for _ in cols], cols)
    for u in keys(Sp)
        b = B[unitchannels[u]]
        is = AN.interval.(b)
        whole = length(Sp[u])
        burst = burst = AN.inany.(Sp[u], (is,)) |> sum
        nonburst = whole - burst

        ΔT = (maximum(Sp[u]) - minimum(Sp[u]))
        Δt = IntervalSets.width.(is) |> sum

        whole = whole/ΔT
        burst = burst/Δt
        nonburst = nonburst/(ΔT - Δt)

        push!(stats, [u, unitchannels[u], whole, burst, nonburst])
    end
    return stats
end
