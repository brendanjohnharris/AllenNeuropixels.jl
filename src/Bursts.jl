using StatsBase
using ImageMorphology
using LsqFit
using LinearAlgebra
using ProgressLogging
using Interpolations
using HypothesisTests

abstract type AbstractBurst end
BurstVector = AbstractVector{<:AbstractBurst}

Base.@kwdef mutable struct Burst <: AbstractBurst
    mask::LogWaveletMatrix
    thresh::Any
    peak = nothing
    fit = nothing
    width = nothing
    significance = nothing
    phasemask::Union{Nothing, LogWaveletMatrix} = nothing
end
Burst(mask, thresh; kwargs...) = Burst(; mask, thresh, kwargs...)
Burst(mask, thresh, peak; kwargs...) = Burst(; mask, thresh, peak, kwargs...)

phasemask(B::Burst) = B.phasemask
duration(B::Burst) = B.fit.param[4]
logspectralwidth(B::Burst) = B.fit.param[5]
maxperiod(B::Burst) = 1 / maxfreq(B)
function periods(B::Burst)
    duration(B) / maxperiod(B)
end
function spectralwidth(B::Burst) # Std of a log-normal distribution
    σ = logspectralwidth(B)
    μ = logpeakfreq(B)
    return sqrt((exp10(σ^2) - 1) * exp10(2 * μ + σ^2))
end
amplitude(B::Burst) = B.fit.param[1]
mask(B::Burst) = B.mask
# binarymask(B::Burst) = B.mask .> B.thresh[1]
peaktime(B::Burst) = B.fit.param[2]
maxmaskfreq(B::Burst) = exp10(dims(mask(B), 2)[end])
minmaskfreq(B::Burst) = exp10(dims(mask(B), 2)[1])
maxmasktime(B::Burst) = dims(mask(B), 1)[end]
minmasktime(B::Burst) = dims(mask(B), 1)[1]
maskfrequencyinterval(B::Burst) = minmaskfreq(B) .. maxmaskfreq(B)
masklogfrequencyinterval(B::Burst) = log10(minmaskfreq(B)) .. log10(maxmaskfreq(B))
maskinterval(B::Burst) = minmasktime(B) .. maxmasktime(B)
maskfreqextrema(B::Burst) = (minmaskfreq(B), maxmaskfreq(B))
maxfreq(B::Burst) = exp10(dims(mask(B), 2)[findmax(mask(B))[2][2]])
maxtime(B::Burst) = dims(mask(B), 1)[findmax(mask(B))[2][1]]
logpeakfreq(B::Burst) = B.fit.param[3]
peakfreq(B::Burst) = B |> logpeakfreq |> exp10
df(B::Burst) = mean(diff(collect(dims(mask(B), Log𝑓))))
dt(B::Burst) = step(dims(mask(B), 𝑡))
maskduration(B::Burst) = dims(mask(B), 𝑡) |> extrema |> collect |> diff |> first
function maskspectralwidth(B::Burst)
    dims(mask(B), Log𝑓) |> extrema |> collect |> diff |> first
end
asymmetry(B::Burst) = abs(log10(maxfreq(B)) - log10(peakfreq(B)))
fiterror(B::Burst) = mean(B.fit.resid) ./ mean(B.mask[:])
interval(B::Burst, σ = 0.5) = peaktime(B) ± (σ * duration(B))
timeinterval = interval
logfrequencyinterval(B::Burst, σ = 0.5) = logpeakfreq(B) ± (σ * logspectralwidth(B))
inany(x::Number, V::Vector{<:AbstractInterval}) = any(in.((x,), V))
inany(x::Vector, V::Vector{<:AbstractInterval}) = inany.(x, (V,))
inany(x::Tuple, V::Vector{<:AbstractInterval}) = inany.(x, (V,))
isfit(B::Burst) = (getfield(B, :fit) isa LsqFit.LsqFitResult)

function getchannel(B::Burst)
    isempty(mask(B).refdims) ? nothing : refdims(mask(B), Chan) |> first
end
flds = [:stimulus, :structure, :sessionid, :probeid]
for f in flds
    fs = "get" * string(f) |> Symbol
    @eval ($fs)(B::Burst) = haskey(B.mask.metadata, Symbol($f)) ?
                            B.mask.metadata[Symbol($f)] : nothing
end
# stimulus(B::Burst) = hasfield(B.mask.metadata, :stimulus) ? B.mask.metadata[:stimulus] : nothing
# structure(B::Burst) = hasfield(B.mask.metadata, :structure) ? B.mask.metadata[:structure] : nothing
# sessionid(B::Burst) = hasfield(B.mask.metadata, :sessionid) ? B.mask.metadata[:sessionid] : nothing
# channel(B::Burst) = hasfield(B.mask.metadata, nnnel) ? B.mask.metadata[:channel] : nothing
# probeid(B::Burst) = hasfield(B.mask.metadata, :probeid) ? B.mask.metadata[:probeid] : nothing
function _burstsubset(B, σ)
    channel = getchannel(B)
    isnothing(channel) && @error "Burst has no associated LFP"
    ts = interval(B, σ)
    return channel, ts
end

function burstsubset(B::Burst, LFP::LFPMatrix, σ = 1.0)
    channel, ts = _burstsubset(B, σ)
    return LFP[𝑡(ts), Chan(At(channel))]
end

function burstsubset(B::Burst, res::LogWaveletMatrix, σ = 1.0)
    channel, ts = _burstsubset(B, σ)
    @assert channel == getchannel(B)
    return res[𝑡(ts)]
end

burstlfp = burstsubset

function basicfilter!(B::BurstVector; pass = nothing, fmin = 0.1, tmin = 2, tmax = 5) # fmin in octaves
    # fmin = log10(fmin)
    # pass = Interval(pass...)
    # passes(x) = all([y in pass for y in x])
    # filter!(b->size(mask(b), 𝑡) > tmin/dt(b), B)
    # filter!(b->size(mask(b), 𝑡) < tmax/dt(b), B)
    filter!(b -> size(mask(b), 𝑡) > tmin / maxfreq(b) / dt(b), B) # Filter bursts that have a duration less than tmin cycles
    isnothing(pass) || (filter!(b -> size(mask(b), 𝑡) < tmax / pass[1] / dt(b), B))
    # filter!(b->size(mask(b), Log𝑓) > fmin/df(b), B)
    # filter!(b->passes(extrema(dims(mask(b), 𝑓))), B)
end
basicfilter!(; kwargs...) = x -> basicfilter!(x; kwargs...)

function bandfilter!(B::BurstVector; pass = [50, 60])
    pass = Interval(pass...)
    filter!(b -> maxfreq(b) ∈ pass, B)
    filter!(b -> peakfreq(b) ∈ pass, B)
end
function periodfilter!(B::BurstVector; n = 1)
    filter!(b -> periods(b) ≥ n, B)
end
function periodfilter!(B::BurstVector, B_sur::BurstVector; σ = 2)
    ps = periods.(B_sur)
    f(b) = periods(b) > median(ps) + iqr(ps) * σ
    filter!(f, B)
end
filterbursts! = bandfilter!

function filtergof!(B::BurstVector)
    # filter!(b -> minmaskfreq(b) < peakfreq(b) < maxmaskfreq(b), B)
    filter!(b -> minmasktime(b) < peaktime(b) < maxmasktime(b), B)
    filter!(b -> maxmasktime(b) - minmasktime(b) > duration(b), B)
    # filter!(b -> maxmaskfreq(b) - minmaskfreq(b) > spectralwidth(b), B)
end

bandfilter!(; kwargs...) = x -> bandfilter!(x; kwargs...)
bandfilter(B; kwargs...) = (Bs = deepcopy(B); bandfilter!(Bs; kwargs...); Bs)
bandfilter(; kwargs...) = x -> bandfilter(x; kwargs...)

function filtergammabursts!(B::BurstVector; fmin = 1, tmin = 0.008, pass = [50, 60]) # fmin in Hz, tmin in s
    fmin = log10(fmin)
    pass = Interval(pass...)
    filter!(b -> size(mask(b), 𝑡) > tmin / dt(b), B)
    filter!(b -> size(mask(b), Log𝑓) > fmin / df(b), B)
    filter!(b -> peakfreq(b) ∈ pass, B)
    # ..... Other stuff
end

filternbg! = filtergammabursts!

# # Old method, global distribution
# function threshold(res, thresh, method)
#     if method == :percentile
#         cutoff = percentile(res[:], thresh)
#     elseif method == :std
#         cutoff = thresh*std(res[:])
#     end
# end
function threshold(res, thresh, method)
    if !(method isa Symbol) # A tuple of (method, surrogate_res)
        sres = last(method)
        sres = convert(LogWaveletMatrix, sres)
        dt = step(dims(res, 𝑡))
        tilims = collect(extrema(dims(res)[1]))
        tilims[2] += dt / 4 # Slightly widen to avoid floting point issues with the indices
        tilims[1] -= dt / 4 # Slightly widen to avoid floting point issues with the indices
        res = sres[𝑡(ClosedInterval(tilims...)),
                   Log𝑓(ClosedInterval(extrema(dims(res)[2])...))]
        method = first(method)
    end
    if method == :percentile
        cutoff = mapslices(x -> percentile(x, thresh), res, dims = 𝑡)
    elseif method == :std
        cutoff = thresh * mapslices(std, res, dims = 𝑡) # .+ mapslices(mean, res, dims = 𝑡)
    elseif method == :iqr
        cutoff = thresh * mapslices(iqr, res, dims = 𝑡) #.+
        mapslices(median, res, dims = 𝑡)
    end
    @debug "Thresholding at $thresh with $(method) method"
    return cutoff
end

# function surrogatewavelettransform(x::LFPVector; method=AP(samplingrate(x)), transform=AN.mmapwavelettransform, kwargs...)
#     y = surrogate(x, method)
#     sres = transform(y; kwargs...)
# end

# grad(x::Real, y::Real) = x - y
# grad(x::AbstractVector, h) = grad.((@view x[3:end]), (@view x[1:end-2]))./(2*h)
# grad(h::Real) = x -> grad(x, h)

"""
Threshold a wavelet spectrum using either a percentile cutoff (`method=:percentile`) or a standard deviation cutoff (`method=:std`) of either each frequency band (`eachfreq=true`) or the entire spectrum. You probably want to FOOOF the spectrum before this
"""
function burstthreshold!(res::LogWaveletMatrix, thresh; method = :iqr, zerograd = 0.0)
    @assert dims(res, 𝑡).val.data isa AbstractRange "Rectify the LFP array before calculating the wavelet transform"
    cutoffs = threshold(res, thresh, method)
    res[res .< cutoffs] .= 0.0
end

function burstthreshold(res, thresh; kwargs...)
    (y = deepcopy(res); burstthreshold!(y, thresh; kwargs...); y)
end

function burstthreshold!(res::WaveletMatrix, thresh; kwargs...)
    burstthreshold!(convert(LogWaveletMatrix, res), thresh; kwargs...)
end

"""

`thresh` is a proportion of the average gradient of the wavelet spectrum below which a gradie
"""
function burstcurvature!(res::LogWaveletMatrix, thresh = 0)
    @assert dims(res, 𝑡).val.data isa AbstractRange "Rectify the LFP array before calculating the wavelet transform"

    xs = dims(res, 𝑡)
    @assert xs.val.data isa AbstractRange "Expected rectified time indices"
    xs = xs.val.data
    ys = collect(dims(res, Log𝑓))
    @assert std(diff(ys)) / std(ys)<1e-3 "Logarithmic frequency bins are not approximately uniform"
    df = median(diff(ys))
    _ys = minimum(ys):df:(maximum(ys) + 3.0 * df)
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
    res[kx .> -thresh .* median(kx)] .= 0.0
    res[ky .> -thresh .* median(ky)] .= 0.0
end
function burstcurvature(res, thresh = 0; kwargs...)
    (y = deepcopy(res); burstcurvature!(y, thresh; kwargs...); y)
end
function burstcurvature!(res::WaveletMatrix, thresh = 0; kwargs...)
    burstcurvature!(convert(LogWaveletMatrix, res), thresh; kwargs...)
end

function widen(x, δ = 0.5; upperbound = [Inf, Inf])
    δ = δ / 2
    @assert length(x[1]) == length(x[2]) == 2
    Δ = [x[2][1] - x[1][1], x[2][2] - x[1][2]]
    return [
            max.(1, floor.(Int, x[1] .- δ .* Δ)),
            min.(upperbound, ceil.(Int, x[2] .+ δ .* Δ))
            ]
end

function _detectbursts(res::LogWaveletMatrix; thresh = 4, curvaturethresh = 2,
                       boundingstretch = 0.5, method = :iqr, areacutoff = 1, dofit = false,
                       filter = nothing)
    @debug "Thresholding amplitudes"
    _res = burstthreshold(res, thresh; method) .> 0
    @debug "Thresholding curvatures"
    _res = _res .& (burstcurvature(res, curvaturethresh) .> 0)

    @debug "Finding connected components"
    components = ImageMorphology.label_components(_res |> Array)
    areas = ImageMorphology.component_lengths(components)[2:end]
    idxs = areas .≥ areacutoff
    centroids = ImageMorphology.component_centroids(components)[2:end]
    centroids = [round.(Int, c) for c in centroids]
    frq = dims(_res, Log𝑓)
    t = dims(_res, 𝑡)
    peaks = [(t[c[1]], frq[c[2]]) for c in centroids]
    peaks = [(p[1], p[2], _res[𝑡(At(p[1])), Log𝑓(At(p[2]))]) for p in peaks]
    # ! Might want to look for gaussian peak instead
    # masks = burstmask.((_res,), peaks)
    bb = ImageMorphology.component_boxes(components)[2:end]
    bb = [widen(b, boundingstretch; upperbound = size(res)) for b in bb]
    masks = [res[b[1][1]:b[2][1], b[1][2]:b[2][2]] for b in bb]
    B = Burst.(masks, ((method, thresh, curvaturethresh),), peaks)[idxs]
    isnothing(filter) || filter(B)

    B = [b for b in B if all(size(b.mask) .> 1)]

    if dofit
        @info "Fitting $(length(B)) bursts"
        fit!(B)
        B = B[isfit.(B)]
    end
    return B
end

function mmap_detectbursts(res::LogWaveletMatrix; window = 50000, kwargs...)
    ti = _slidingwindow(collect(dims(res, 𝑡)), window; tail = true)
    ti = [ClosedInterval(extrema(t)...) for t in ti]
    B = Vector{Burst}()
    threadlog, threadmax = (0, length(ti))
    @withprogress name="Burst detection" begin
        for (i, ts) in enumerate(ti)
            @info "Calculating epoch $i/$(length(ti))"
            subres = res[𝑡(ts)]
            append!(B, _detectbursts(subres; kwargs...))
            if threadmax > 1
                Threads.threadid() == 1 && (threadlog += 1) % 1 == 0 &&
                    @logprogress threadlog / threadmax
            end
        end
    end
    return B
end

function pvalues(B::BurstVector, Bs::BurstVector, f::Function; test = OneSampleTTest,
                 tail = :left)
    d = f.(B)
    ds = f.(Bs)
    p = HypothesisTests.pvalue.(test.((ds,), d); tail) .* length(d) # Bonferroni correction
end

function pvalues(B::BurstVector, Bs::BurstVector, fs::AbstractVector{<:Function};
                 test = OneSampleHotellingT2Test, tail = nothing)
    d = collect.(zip([f.(B) for f in fs]...))
    ds = hcat([f.(Bs) for f in fs]...)
    tests = test.((ds,), d)
    p = HypothesisTests.pvalue.(tests) .* length(d) # Bonferroni correction
end

function significancefilter!(B::BurstVector, Bs::BurstVector, f; thresh = 0.05, kwargs...)
    p = pvalues(B, Bs, f; kwargs...)
    deleteat!(B, p .> thresh)
end

function mahalanobisfilter!(B::BurstVector, Bs::BurstVector; σ = 2,
                            metrics = [amplitude, periods])
    s = hcat([m.(Bs) for m in metrics]...)
    Σ² = cov(s)
    μ = mean(s, dims = 1)'
    x = hcat([m.(B) for m in metrics]...)
    ℳ(x) = sqrt((x - μ)' * inv(Σ²) * (x - μ))
    p = ℳ.(eachrow(x)) .|> first
    # idxs = [_x .> mean(_s) for (_x, _s) in zip(eachcol(x), eachcol(s))] # Ensure we only take the ones above the respective means
    # idxs = idxs[1] .| idxs
    deleteat!(B, p .< σ)
end
function bivariatefilter!(B::BurstVector, Bs::BurstVector; 𝑝 = 0.05,
                          metrics = [amplitude, periods])
    s = hcat([m.(Bs) for m in metrics]...)
    f(p) = HypothesisTests.pvalue(OneSampleHotellingT2Test(s, [m(p) for m in metrics])) < 𝑝
    filter!(f, B)
end

function metricfilter!(B::BurstVector, Bs::BurstVector, metric; σ = 3)
    s = metric.(Bs)
    thresh = median(s) .+ iqr(s) * σ
    filter!(b -> metric(b) > thresh, B)
end

"""
Filter a vector of bursts based of a vector of surrogate bursts
"""
function surrogatefilter!(B::BurstVector, Bs::BurstVector; stats = [duration, peakfreq],
                          joint = true)
    if joint
        significancefilter!(B, Bs, stats; test = OneSampleHotellingT2Test, tail = nothing)
    else
        for s in stats
            # ! The problem is, the distribution changes with each of these. Will want to take the maximum p-value, rather successively filtering
            significancefilter!(B, Bs, s; test = OneSampleTTest, tail = :left,
                                thresh = 0.01)
        end
    end
end

"""
Detect bursts from a supplied wavelet spectrum, using thresholding
`boundingstretch` increases the bounding box slightly so for a more accurate fit. Give as a proportion of the threshold bounding box
`detection` can be `_detectbursts` or `mmap_detectbursts`
"""
function detectbursts(res::LogWaveletMatrix; pass = nothing, dofit = true,
                      detection = _detectbursts, kwargs...)
    if !isnothing(pass)
        fpass = collect(log10.(pass))
        fpass = [fpass[1] - 0.25, fpass[2] + 0.25] # Slightly widen for better burst estimates. We still filter burst centers in the provided pass band later
        @info "Selecting the $(pass) Hz band"
        res = res[Log𝑓(Interval(fpass...))]
    end
    B = detection(res; filter = basicfilter!(; pass), dofit, kwargs...)
    # isnothing(pass) || (@info "Filtering in the $(pass) Hz band"; bandfilter!(B; pass))

    if dofit
        sort!(B, by = peaktime)
    end
    return B
end

function detectbursts(x::LFPVector; kwargs...) # surrodur=min(length(x), round(Int, 50/step(dims(x, 𝑡))/minimum(pass))), N=100
    # s = surrogate(x, AP())
    # # S = [s() for _ in 1:N]
    # # γₛ = gammafilter(s; pass)
    # res = wavelettransform(γₛ)
    # agg = [getindex.(res, (i,)) for i in CartesianIndices(res[1])]
    # agg = ToolsArray(agg, dims(res[1]))
    # σ = mapslices(std, res, dims=1)
    # μ = mapslices(mean, res, dims=1)

    # γ = gammafilter(x; pass)
    res = fooofedwavelet(x)
    detectbursts(res; kwargs...)
end

function simpleburstdetection(y_back; thresh = 5, pass, method = :iqr, kwargs...)
    # First construct the group estimate and the surrogate
    res_back = _wavelettransform(y_back)
    # res_back = AN.logwaveletmatrix(res_back .|> abs)
    phi = angle.(res_back)
    res = abs.(res_back)
    res = logwaveletmatrix(res)
    phi = logwaveletmatrix(phi)
    sres = surrogate(y_back, FT()) |> _wavelettransform |> logwaveletmatrix
    sphi = angle.(sres)
    sres = abs.(sres)
    B = detectbursts(res; detection = detectbursts, dofit = true,
                     pass, method, thresh, kwargs...)
    addphasemasks!(B, phi)
    bandfilter!(B; pass)
    filtergof!(B)

    B_sur = detectbursts(sres; detection = detectbursts, dofit = true,
                         pass, method, thresh, kwargs...)
    addphasemasks!(B_sur, sphi)
    bandfilter!(B_sur; pass)
    filtergof!(B_sur)

    @assert length(B_sur)>5 "Not enough surrogate bursts to estimate a distribution"
    mahalanobisfilter!(B, B_sur; σ = 2, metrics = [amplitude, periods])
    # metricfilter!(B, B_sur, amplitude; σ = 1)

    periodfilter!(B, n = 1)
    periodfilter!(B_sur, n = 1)
    return B, B_sur
end

function fit!(B::Burst)
    mask = B.mask

    # The fitting takes a while. Let's downsample a little bit if we can
    while size(mask, 1) > 80
        mask = mask[1:2:end, :]
    end
    while size(mask, 2) > 80
        mask = mask[:, 1:2:end]
    end
    B.fit = fitdiagonalgaussian(mask)
end

function fit!(ℬ::AbstractVector{<:AbstractBurst})
    @withprogress name="Fitting bursts" begin
        threadlog, threadmax = (0, length(ℬ))
        l = Threads.ReentrantLock()
        Threads.@threads for i in eachindex(ℬ)
            try
                fit!(ℬ[i])
            catch e
                @warn e
                @warn "Could not fit a burst. Skipping"
            end
            lock(l)
            try
                threadlog += 1
                if threadlog > threadmax ÷ 3
                    @logprogress threadlog / threadmax
                end
                # @debug "Burst fitting progress: $threadlog/$threadmax"
            finally
                unlock(l)
            end
        end
    end
end

function gaussian2!(F, xy, p)
    x, y = eachcol(xy)
    A, μ₁, μ₂, σ₁, σ₂ = p
    @. F = A * exp(-0.5 * (((x - μ₁) / σ₁)^2 + ((y - μ₂) / σ₂)^2))
end
gaussian2(xy::Tuple, p) = (F = [0.0]; gaussian2!(F, collect(xy)', p) |> first)

function gaussian2_j!(J::Array{Float64, 2}, xy, p)
    x, y = eachcol(xy)
    A, μ₁, μ₂, σ₁, σ₂ = p
    @. J[:, 1] .= exp(-0.5 * (((x - μ₁) / σ₁)^2 + ((y - μ₂) / σ₂)^2))                       #dF/A
    @. @views J[:, 2] .= A * (x - μ₁) * J[:, 1] / (σ₁^2)         #dF/μ₁
    @. @views J[:, 3] .= A * (y - μ₂) * J[:, 1] / (σ₂^2)         #dF/μ₂
    @. @views J[:, 4] .= A * (x - μ₁)^2 * J[:, 1] / (σ₁^3)       #dF/σ₁
    @. @views J[:, 5] .= A * (x - μ₂)^2 * J[:, 1] / (σ₂^3)       #dF/σ₂
end

function fitdiagonalgaussian(x, y, Z::AbstractMatrix)
    x, y = collect.((x, y))
    xy = hcat(collect.(Iterators.product(x, y))...)'
    z = Z[:]
    A0, μ0 = findmax(z)
    p0 = [A0, # A
        xy[μ0, 1], # μ₁
        xy[μ0, 2], # μ₂
        (x |> extrema |> collect |> diff |> first) / 3, # σ₁
        (y |> extrema |> collect |> diff |> first) / 3] # σ₂
    T = eltype(z)
    xy = xy .|> T
    p0 = p0 .|> T
    fit = LsqFit.curve_fit(gaussian2!, xy, z, p0; inplace = true, autodiff = :forwarddiff)
end

function fitdiagonalgaussian(mask::AbstractToolsArray)
    frq = dims(mask, 2)
    t = dims(mask, 1)
    f = fitdiagonalgaussian(t, frq, mask |> Array)
end

# """
# Get a sufficient image of a burst from an entire wavelet spectrum. A threshold at which the power is 1% of the peak power should suffice, but really the threshold should just be small enough for the exact value not to matter.
# Inputs are `res`, the complete wavelet spectrum, `peak`, the (time, freq, value) or (time, freq) tuple of the peak, and `thresh`, the proportion threshold on the peak power. `diffthresh` is a threshold on the gradient of the spectrum, to avoid capturing a bimodal burst.

# Could consider doing this in two passess, to be computationally efficient: first pass by looking along perpendicular directions, with a low threshold, then second pass by fitting a Gaussian on this loose mask and setting a higher threshold.
# """
function burstmask(res::LogWaveletMatrix, peak; thresh = 0.8, n = 3)#, diffthresh=0.01)
    pidx = ((𝑡 ∘ At)(peak[1]), (Log𝑓 ∘ At)(peak[2]))
    A = res[pidx...]
    thresh = thresh * A
    # Find the rectangular boundary that contains the threshold. Work in each dimension from the centre. Assumes the burst profile is blobby around the peak and its horizontal/vertical axes.
    frq = dims(res, Log𝑓)
    t = dims(res, 𝑡)

    bounds = zeros(4) # (tmin, tmax, fmin, fmax)

    # ! This method is not that great
    # ups = res[pidx[2]][t .< peak[1]]
    # dups = -vcat(diff(ups), [0])
    # up = findlast((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[1] = isnothing(up) ? 1 : up

    # ups = res[pidx[2]][t .> peak[1]]
    # dups = vcat([0], diff(ups))
    # up = findfirst((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[2] = isnothing(up) ? size(res, 𝑡) : findfirst(t .> peak[1]) + up - 1

    # ups = res[pidx[1]][frq .< peak[2]]
    # dups = -vcat(diff(ups), [0])
    # up = findlast((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[3] = isnothing(up) ? 1 : up

    # ups = res[pidx[1]][frq .> peak[2]]
    # dups = vcat([0], diff(ups))
    # up = findfirst((ups .< thresh) .| (dups .> diffthresh*maximum(dups)))
    # bounds[4] = isnothing(up) ? size(res, 𝑓) : findfirst(frq .> peak[2]) + up - 1

    # bounds[1:2] .= t[Int.(bounds[1:2])]
    # bounds[3:4] .= frq[Int.(bounds[3:4])]

    # * This one looks for the hwhm in each direction, then just multiplies that by a few for good coverage. It's not a terrible thing if we are too broad.
    p = t .< peak[1]
    ups = res[pidx[2]][p]
    up = findlast((ups .< thresh))
    w = length(ups) - up
    up = findlast(p) - n * w
    bounds[1] = max(1, up)

    p = t .> peak[1]
    ups = res[pidx[2]][p]
    up = findfirst((ups .< thresh))
    w = up
    up = findfirst(p) + n * w
    bounds[2] = min(size(res, 𝑡), up)

    p = frq .< peak[2]
    ups = res[pidx[1]][p]
    up = findlast((ups .< thresh))
    w = length(ups) - up
    up = findlast(p) - n * w
    bounds[3] = max(1, up)

    p = frq .> peak[2]
    ups = res[pidx[1]][p]
    up = findfirst((ups .< thresh))
    w = up
    up = findfirst(p) + n * w
    bounds[4] = min(size(res, 𝑓), up)

    bounds[1:2] .= t[Int.(bounds[1:2])]
    bounds[3:4] .= frq[Int.(bounds[3:4])]

    mask = res[𝑡(Interval(bounds[1:2]...)), 𝑓(Interval(bounds[3:4]...))]
end

"""
Pair each burst in A with a burst in B, i.e. where the burst in A is closest in time to the burst in B. Not all bursts will be paired (?)
First, calculate a distance matrix between all bursts.
Then, select bursts pairs with a minimum distance from the distance matrix
"""
# function pairbursts(A::BurstVector, B::BurstVector; feed=:forward)
#     B = sort(deepcopy(B), by=peaktime)
#     A = sort(deepcopy(A), by=peaktime)
#     C = []
#     D = []
#     ta = peaktime.(A)
#     tb = peaktime.(B)
#     for (i, a) in enumerate(A)
#         ts = ta[i] .> tb
#         t = findlast(ts)
#         if ~isnothing(t)
#             if feed == :forward
#                 if t+1 > length(ts)
#                     break
#                 end
#                 ts[t+1] = true
#                 t = findlast(ts)
#             end
#             push!(C, a)
#             push!(D, B[t])
#             deleteat!(B, ts) # Remove all bursts before `a``
#             deleteat!(tb, ts)
#         end
#     end
#     return C, D
# end

"""
May have duplicate bursts...
"""
# function pairdistances(Δ::Matrix; side=:both)
#     cols = []
#     for i in axes(Δ, 1)
#         _, d = findmin(abs.(Δ[i, :])) # Find the burst in the secondary region with the smallest distance to the primary burst
#         push!(cols, d)
#     end
#     Δ = Δ[:, cols]
#     if side == :both
#         _cols = []
#         _rows = []
#         for c in unique(cols)
#             _, I = findmin(abs.(Δ[:, cols .== c]))
#             i, j = I |> Tuple
#             push!(_cols, j)
#             push!(_rows, i)
#         end
#         Δ = Δ[_rows, _cols]
#         cols = (_rows, cols)
#     end
#     return Δ, cols
# end

"""
Detect paired bursts that are closer than thresh
"""
# function pairdistances(Δ, thresh=0.25)
#     I = findall(abs.(Δ) .< thresh)
#     I = I .|> Tuple
#     i = unique(first.(I))
#     j = unique(last.(I))
#     return Δ[i, j], (i, j)
# end
# function pairdistances(Δ)
#     Δ = deepcopy(Δ)
#     # * Step along rows, finding the closest pairs of bursts
#     _Δ = fill(NaN, (size(Δ, 1), size(Δ, 1)))
#     while any(isnan.(_Δ))
#         r, c = Δ .|> abs |> findmin |> last |> Tuple
#         _Δ[:, r] .= Δ[:, c]
#         Δ[:, c] .= Inf
#         Δ[r, :] .= Inf
#     end
#     return _Δ
# end
function pairdistances(Δ)
    I = [abs.(x) |> findmin |> last for x in eachrow(Δ)]
    return Δ[:, I] # The bursts in reg. 1, and the distance to their closes bursts in reg. 2
end

function burstdistances(A::BurstVector, B::BurstVector)
    ta = peaktime.(A)
    tb = peaktime.(B)
    return [b - a for a in ta, b in tb]
end

function burstdelta(A::BurstVector, B::BurstVector; feed = :forward)
    # C, D = pairbursts(A, B; feed)
    # Δ = peaktime.(D) .- peaktime.(C)
    Δ = burstdistances(A, B)
    I = [x |> findmin |> last for x in abs.(Δ) |> eachrow]
    return [x[i] for (x, i) in zip(eachrow(Δ), I)]
end

function burstspikestats(B, Sp, channels; sessionid, probeid, phifreqs = 1:1:100, kwargs...)
    phifreqs = log10.(phifreqs)
    unitchannels = getclosestchannels(sessionid, probeid, keys(Sp), channels)
    # * Count the number of spikes within bursts vs outside bursts. Assume Sp contains spikes only in the duration of the LFP used to calculate bursts.
    cols = [
        :unit,
        :channel,
        :mean_rate,
        :burst_rate,
        :nonburst_rate,
        :rate_index,
        :phase_synchrony
    ]
    # if !isnothing(phi)
    #     append!(cols, [:phase_synchrony])
    # end
    stats = DataFrame([[] for _ in cols], cols)
    for u in keys(unitchannels)
        b = B[unitchannels[u]]
        is = interval.(b)
        whole = length(Sp[u])
        if isempty(is)
            topush = [u, unitchannels[u], fill(NaN, length(cols) - 2)...]
        else
            _burst = inany.(Sp[u], (is,))
            burstspikes = Sp[u][_burst]
            burst = _burst |> sum
            nonburst = whole - burst

            ΔT = (maximum(Sp[u]) - minimum(Sp[u]))
            Δt = IntervalSets.width.(is) |> sum

            whole = whole / ΔT
            burst = burst / Δt
            nonburst = nonburst / (ΔT - Δt)

            rate_index = burst ./ nonburst

            phase_synchrony = phaselockingindex.([b], [Sp[u]], phifreqs)
            # phaselockingindex(phi, burstspikes) # burstspikes are only the spikes that occur during bursts
            topush = [
                Int(u),
                unitchannels[u],
                whole,
                burst,
                nonburst,
                rate_index,
                phase_synchrony
            ]
        end
        push!(stats, topush)
    end
    return stats
end

function _phaselockingindex(phi::LogWaveletMatrix, s::AbstractVector)
    # * Check the spikes are all in the bounds
    inter = ClosedInterval(extrema(dims(phi, 𝑡))...)
    s = s[s .∈ [inter]]
    # * Calculate the phases during each spike, for each frequency
    phis = phi[𝑡(Near(s))]
    return phis
end

function pairwisephaseconsistency(x::AbstractVector) # Eq. 14 of Vinck 2010
    # x = collect(x)
    N = length(x)
    # f(ϕ, ω) = cos(ϕ - ω) # cos(ϕ)*cos(ω) + sin(ϕ)*sin(ω) # Dot product between unit vectors with given phases
    # f(a) = f(a...)
    Δ = zeros(N - 1)
    Threads.@threads for i in 1:(N - 1)
        δ = @views x[i] .- x[(i + 1):end]
        Δ[i] = sum(cos.(δ))
    end
    # Δ = (sum(f.(b)) - N)/2 # -N to remove the diagonal of 1's, /2 to remove lower triangle
    return (2 / (N * (N - 1))) * sum(Δ)
end
function pairwisephaseconsistency(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y)
    z = y .- x
    return pairwisephaseconsistency(z)
end

"""
Buzsaki's phase-locking index ("Gamma rhythm communication between entorhinal cortex and dentate gyrus neuronal assemblies")
"""
function phaselockingindex(phi::LogWaveletMatrix, s::AbstractVector)
    phis = _phaselockingindex(phi, s)
    𝒴 = mapslices(pairwisephaseconsistency, phis, dims = 𝑡)[𝑡(1)]
end

"""
Calculate the phase-locking index using the wavelet transform masks stored in each bursts. Drops any bursts that do not have wavelet information at the specified frequency, `f`. If `centre` is true, the phase of spikes in each burst is centred to the mean phase in each burst. Only consider bursts with at least `n_spikes` spikes, and only count channel-neuron pairs that have `n_bursts` bursts.
"""
function _phaselockingindex(B::BurstVector, s::AbstractVector, f::Number; centre = true,
                            n_spikes = 5, n_bursts = 10)
    # First, get a list of phases for every spike with the burst interval
    ts = interval.(B)
    phis = phasemask.(B)
    phis = getindex.(phis, [Log𝑓(Near(log10(f)))])
    s = s[inany(s, ts)]
    ϕ = [zeros(0) for _ in phis]
    for (e, es) in enumerate(s)
        i = findfirst([es ∈ t for t in ts])
        ϕ_ = phis[i][𝑡(Near(es))]
        append!(ϕ[i], ϕ_)
    end
    if centre
        ϕ = [p .- (im .* p .|> exp |> sum |> angle) for p in ϕ] # Subtract circular mean
    end
    if n_spikes > 0
        ϕ = [p for p in ϕ if length(p) ≥ n_spikes]
    end
    N_bursts = length(ϕ)
    if n_bursts > 0 && N_bursts < n_bursts
        # @warn "This channel has an insufficient number of bursts with spikes: $(length(ϕ))"
        return []
    end
    ϕ = vcat(ϕ...)
    # @info "This channel-spike pair has $(N_bursts) bursts and $(length(ϕ)) spikes"
    return ϕ
end

"""
The phase-locking index for bursts and an LFP signal
"""
function _phaselockingindex(B::BurstVector, s::LFPVector, f::Number; centre = true,
                            n_bursts = 5)
    # First, get phases for every burst
    ts = interval.(B)
    phis = phasemask.(B)
    phis = getindex.(phis, [Log𝑓(Near(log10(f)))])

    phis = [p[𝑡(t)] for (p, t) in zip(phis, ts)] # Just to be sure, same span as below
    ts = [dims(p, 𝑡) |> collect for p in phis] # Just to be consistent
    phis′ = [s[𝑡(Near(t))] for t in ts]
    @assert length.(phis) == length.(phis′)

    # Now get phase differences within each burst
    ϕ = [collect(p .- p′) for (p, p′) in zip(phis, phis′)]

    if centre
        ϕ = [p .- (im .* p .|> exp |> sum |> angle) for p in ϕ] # Subtract circular mean
    end
    N_bursts = length(ϕ)
    if n_bursts > 0 && N_bursts < n_bursts
        @debug "This channel has an insufficient number of bursts with spikes: $(length(ϕ))"
        return []
    end
    ϕ = vcat(ϕ...)
    @debug "This channel-spike pair has $(N_bursts) bursts and $(length(ϕ)) spikes"
    return ϕ
end

function phaselockingindex(B::BurstVector, s::AbstractVector, f::Number)
    phis = _phaselockingindex(B, s, f)
    γ = pairwisephaseconsistency(phis)
    𝑝 = isempty(phis) ? 1.0 : HypothesisTests.pvalue(RayleighTest(phis))
    return (γ, 𝑝)
end

function phaselockingindex(LFP::LFPVector, s::AbstractVector; kwargs...)
    phis = LFP |> hilbert .|> angle
    phis = phis[𝑡(Near(s))]
    γ = pairwisephaseconsistency(phis)
    𝑝 = isempty(phis) ? 1.0 : HypothesisTests.pvalue(RayleighTest(collect(phis)))
    return (γ, 𝑝)
end

function phaselockingindex(B::BurstVector, s::LFPVector, f::Number; kwargs...)
    phis = collect(_phaselockingindex(B, s, f; kwargs...))
    γ = pairwisephaseconsistency(phis)
    𝑝 = isempty(phis) ? 1.0 : HypothesisTests.pvalue(RayleighTest(phis))
    return (γ, 𝑝)
end

function phaselockingindex(ℬ::Dict, Sp::Dict, f::Number; kwargs...)
    channels = keys(ℬ) |> collect
    units = keys(Sp) |> collect
    γ = ToolsArray(collect(zeros(length(channels), length(units))),
                   (Chan(channels), Unit(units)))
    𝑝 = deepcopy(γ)
    Threads.@threads for (i, b) in collect(enumerate(values(ℬ)))
        for (j, s) in enumerate(values(Sp))
            # @info "Calculating ($i, $j) of $(size(γ))"
            _γ, _𝑝 = phaselockingindex(b, s, f; kwargs...)
            γ[i, j] = _γ
            𝑝[i, j] = _𝑝
        end
    end
    return γ, 𝑝
end

function phaselockingindex(ℬ::Dict, phi::Dict{T, LogWaveletMatrix} where {T}, f::Number;
                           kwargs...)
    channels = keys(ℬ) |> collect
    units = keys(phi) |> collect
    γ = ToolsArray(collect(zeros(length(channels), length(units))),
                   (Chan(channels), Chan(units)))
    𝑝 = deepcopy(γ)
    @withprogress name="LFP-LFP phase-locking index" begin
        threadlog, threadmax = (0, length(values(ℬ)))
        Threads.@threads for (i, b) in collect(enumerate(values(ℬ)))
            for (j, s) in enumerate(values(phi))
                # @info "Calculating ($i, $j) of $(size(γ))"
                _γ, _𝑝 = phaselockingindex(b, s[:, Near(log(f))], f; kwargs...)
                γ[i, j] = _γ
                𝑝[i, j] = _𝑝
            end
            if threadmax > 1
                Threads.threadid() == 1 && (threadlog += Threads.nthreads()) % 10 == 0 &&
                    @logprogress threadlog / threadmax
            end
        end
    end
    return γ, 𝑝
end

function _phaselockingindex(x::LFPVector, y::LFPVector; pass = nothing)
    isnothing(pass) || (x, y = bandpass.((x, y); pass))
    ϕ_x = x |> hilbert .|> angle
    ϕ_y = y |> hilbert .|> angle
    return ϕ_y .- ϕ_x
end

function phaselockingindex(x::LFPVector, y::LFPVector; pass = nothing)
    phis = collect(_phaselockingindex(x, y; pass))
    γ = pairwisephaseconsistency(phis)
    𝑝 = isempty(phis) ? 1.0 : HypothesisTests.pvalue(RayleighTest(phis))
    return γ, 𝑝
end

function phaselockingindex(x::LFPVector, y::LFPVector, N::Number; pass = nothing)
    @assert length(x) == length(y)
    idxs = 1:N:length(x)
    X = [@view(x[i:(i + N - 1)]) for i in idxs if i + N - 1 ≤ length(x)]
    Y = [@view(y[i:(i + N - 1)]) for i in idxs if i + N - 1 ≤ length(y)]
    γ = zeros(length(X))
    𝑝 = zeros(length(X))
    for i in eachindex(γ)
        phis = collect(_phaselockingindex(X[i], Y[i]; pass))
        γ[i] = pairwisephaseconsistency(phis)
        𝑝[i] = isempty(phis) ? 1.0 : HypothesisTests.pvalue(RayleighTest(phis))
    end
    return γ, 𝑝
end

# * We want to detect durations where the theta bursts occur all across the channels
# * To do this, we look at each time point how many channels have an ongoing theta burst. If this is greater than some value, we count this time as a part of the global burst
function globalbursts(B::Dict; thresh = 0.5)
    ts = [extrema.(interval.(b, 1)) for b in values(B)]
    ts = vcat(collect.(vcat(collect.([vcat(e...) for e in ts])...))...)
    # int = Interval(extrema(ts)...)
    ts = minimum(ts):0.005:maximum(ts) # Down sample a little so that we don't take forever
    gbi = globalburstindex(ts, B)
    gbs = (gbi .> thresh) |> ImageMorphology.label_components |>
          ImageMorphology.component_indices
    gbs = gbs[2:end] .|> extrema
    gbs = [Interval(ts[g[1]], ts[g[2]]) for g in gbs]
end

function globalburstindex(t::Number, ints::AbstractVector)
    counts = [inany(t, i) for i in ints]
    return sum(counts) / length(ints)
end

function globaburstindex(ts, B::AbstractVector)
    ints = [interval.(b, 1) for b in B]
    return globalburstindex.(ts, [ints])
end

globalburstindex(ts, B::Dict) = globaburstindex(ts, collect(values(B)))

"""
Returns a function that evaluates the Gaussian fits of a collection of bursts at any time and log-frequency value
"""
function aggregatefit(ℬ::BurstVector; span = 3)
    ts = timeinterval.(ℬ, span) # Relevant bursts are within 2 SD's of the burst centre
    fs = frequencyinterval.(ℬ, span)
    function G(t, f)
        idxs = t .∈ ts
        _ℬ = ℬ[idxs]
        isempty(_ℬ) && return 0.0
        _fs = fs[idxs]
        _fs = [_fs...]
        _ℬ = _ℬ[f .∈ _fs] # Relevant bursts for this time and frequency
        isempty(_ℬ) && return 0.0
        return sum([gaussian2((t, f), b.fit.param) for b in _ℬ]) # Sum of bursts at this value
    end
end

function gaussianmask(B...; ts = nothing, fs = nothing, span = 3)
    ℬ = B[1]
    if isnothing(ts)
        Δt = dt(ℬ[1])
        ts = vcat(collect.(extrema.(timeinterval.(ℬ, span)))...) |> extrema
        ts = ts[1]:Δt:ts[2]
    end
    if isnothing(fs)
        Δf = df(ℬ[1])
        fs = vcat(collect.(extrema.(frequencyinterval.(ℬ, span)))...) |> extrema
        fs = fs[1]:Δf:fs[2]
    end
    masks = []
    for ℬ in B
        mask = ToolsArray(zeros(length(ts), length(fs)), (𝑡(ts), Log𝑓(fs)))
        for b in ℬ
            tt = dims(mask, 𝑡)
            tt = tt[tt .∈ (timeinterval(b, span),)]
            ff = dims(mask, Log𝑓)
            ff = ff[ff .∈ (frequencyinterval(b, span),)]
            g = (x, y) -> gaussian2((x, y), b.fit.param)
            for t in tt
                for f in ff
                    mask[𝑡(At(t)), Log𝑓(At(f))] = g(t, f)
                end
            end
        end
        push!(masks, mask)
    end
    return masks
end

function burstoverlap(res1::LogWaveletMatrix, res2::LogWaveletMatrix; normdims = 𝑓,
                      globalnorm = false)
    normall(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x))
    function normfrequency(x)
        (x .- minimum(x, dims = 𝑡)) ./ (maximum(x; dims = 𝑡) - minimum(x; dims = 𝑡))
    end
    normdims == :all && (res1, res2 = normall.([res1, res2]))
    normdims == 𝑓 && (res1, res2 = normfrequency.([res1, res2]))
    normdims == :frequencymax &&
        (res1, res2 = maximum.([res1, res2], dims = Log𝑓))

    A = globalnorm ? sum(res1, dims = 𝑡) .* sum(res2, dims = 𝑡) : 1.0
    (A .* sum(res1 .* res2, dims = 𝑡))[1, :]
end

burstoverlap(B1::BurstVector, B2::BurstVector) = burstoverlap(gaussianmask(B1, B2)...)

# function burstoverlap(B1::AbstractBurst, B2::AbstractBurst)
#     M1, M2 = gaussianmask([B1], [B2])
#     M1 = M1./sum(M1)
#     M2 = M2./sum(M2)

#     # * Returns the overlap index of the two bursts
# end

function overlapintervals(B1::BurstVector, B2::BurstVector; thresh = 0.1,
                          durationthresh = 0.0) # Threshold on something. Start with the delta between burst centres. 50 ms. Duration of bursts around 100ms
    B1 = B1[duration.(B1) .> durationthresh]
    B2 = B2[duration.(B2) .> durationthresh]
    Δ = abs.(burstdistances(B1, B2))
    # * Select bursts closer than thresh...
    bis = [[], []]
    for bi in 1:size(Δ, 1)
        td = Δ[bi, :]
        if any(td .< thresh)
            push!(bis[1], B1[bi])
            push!(bis[2], B2[findmin(td)[2]])
        end
    end
    return bis
    # ints = []
    # for i in 1:length(bis[1])
    #     b1 = B1[bis[1][i]]
    #     b2 = B2[bis[2][i]]
    #     int1 = interval(b1)
    #     int2 = interval(b2)
    #     int = int1 ∩ int2
    #     !isempty(int) && push!(ints, int)
    # end
    # return ints
end

function burstcommunication(LFPs, ℬ; metric = :weighted_phase_lag_index) # Both args should be 2-vectors.
    # * Find intervals of nearby bursts
    ℬ = overlapintervals(ℬ...)
    ints = [interval.(B) for B in ℬ]
    channels = [getchannel.(B) for B in ℬ]
    ts = [[LFPs[i][𝑡(int), Chan(At(c))] for (c, int) in zip(channels[i], ints[i])]
          for i in 1:length(LFPs)]

    # * Calculate a vector of e.g. PLI's between the two vectors of bursts and their LFP traces
    C = AA.multitaper(ts...)
    W = AA.spectralconnectivity(C; metric)
    return W
end

function burstcommunication(LFPs, ℬ, f; kwargs...)
    W = burstcommunication(LFPs, ℬ; kwargs...)
    Ws = [w[𝑓(Near(Float64(f)))] for w in W]
    return Ws
end

"""
Shuffle times of a burst vector
"""
function randomisebursts(BS::BurstVector)
    ℬ = deepcopy(BS)
    it = extrema(peaktime.(ℬ))
    ts = (rand(length(ℬ)) .* (it[2] - it[1])) .+ it[1] # Random burst times
    # * Adjust the centre parameter
    [B.fit.param[2] = ts[i] for (i, B) in enumerate(ℬ)]
    # * Subtract the centre from the mask dimensions
    for (i, B) in enumerate(ℬ)
        t0 = B.peak[1]
        _t = B.mask.dims[1].val.data .- t0 .+ ts[i]
        B.mask = ToolsArray(B.mask.data, (𝑡(_t), B.mask.dims[2]))
        B.peak = (B.peak[1] - t0 + ts[i], B.peak[2], B.peak[3])
        if hasfield(typeof(B), :phasemask)
            B.phasemask = ToolsArray(B.phasemask.data, (𝑡(_t), B.phasemask.dims[2]))
        elseif ndims(B.significance) == 2 # backards compat
            B.significance = ToolsArray(B.significance.data,
                                        (𝑡(_t), B.significance.dims[2]))
        end
    end
    # [B.mask.dims[1].val.data = B.mask.dims[1].val.data .- ts[i] for (i, B) in enumerate(ℬ)]
    # [B.peak[1] -= ts[i] for (i, B) in enumerate(ℬ)]
    return ℬ
end

function addphasemask!(b::AbstractBurst, ϕ::LogWaveletMatrix)
    m = mask(b)
    tis = dims(m, 𝑡)
    fs = dims(m, Log𝑓)
    b.phasemask = ϕ[𝑡(Near(collect(tis))), Log𝑓(Near(collect(fs)))]
end
addphasemasks!(B::BurstVector, ϕ::LogWaveletMatrix) = addphasemask!.(B, [ϕ])

function addphasemasks(B::BurstVector, ϕ::LogWaveletMatrix)
    B = deepcopy(B)
    for b in B # Just in case these are old bursts, without the phasemask field
        b = Burst(b.mask, b.thresh, b.peak, b.fit, b.width, b.significance)# , nothing
    end
    addphasemasks!(B, ϕ)
    return B
end

function hellinger(A::AbstractBurst, B::AbstractBurst; samechannelallowed = true)
    μ₁ = A.fit.param[2]
    σ₁ = A.fit.param[4]
    μ₂ = B.fit.param[2]
    σ₂ = B.fit.param[4]
    d = abs(μ₁ - μ₂) / 5
    samechannel = false
    if !samechannelallowed
        samechannel = getchannel(A) == getchannel(B)
    end
    if (d > σ₁ && d > σ₂) || samechannel
        H = 1
    else
        s = σ₁^2 + σ₂^2
        H² = 1 - sqrt((2 * σ₁ * σ₂) / (s)) * exp((-1 / 4) * ((μ₁ - μ₂)^2 / s))
        H = sqrt(H²)
    end
    return H
end

function hellinger(A::BurstVector, B::BurstVector)
    H = ones(length(A), length(B))
    Threads.@threads for I in CartesianIndices(H)
        i, j = Tuple(I)
        if i > j # Only fill in the lower triangle
            H[i, j] = hellinger(A[i], B[j])
        elseif i == j
            H[i, j] = 0.0
        end
    end
    H = Symmetric(H, :L)
    return H
end

hellinger(B::BurstVector) = hellinger(B, B)
