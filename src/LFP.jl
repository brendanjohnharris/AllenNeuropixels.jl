using DimensionalData
using IntervalSets
using HDF5
using Statistics
using DSP
using FFTW
using DelayEmbeddings
using MultivariateStats
using Wavelets
using ContinuousWavelets
using TimeseriesSurrogates
using HTTP
using PyFOOOF
using ProgressLogging
using Mmap

bandpass(; kwargs...) = x -> bandpass(x; kwargs...)

function bandpass(x::LFPVector; pass, designmethod = Butterworth(4))
    t = dims(x, Ti)
    T = t isa AbstractRange ? step(t) : t |> collect |> diff |> mean
    fs = 1.0 / T
    y = filtfilt(digitalfilter(Bandpass(pass...; fs), designmethod), x)
    return DimArray(y, dims(x))
end

function bandpass(X::LFPMatrix, dim = Dim{:channel}; kwargs...)
    Y = similar(X)
    Threads.@threads for x in dims(X, dim)
        Y[dim(At(x))] .= bandpass(X[dim(At(x))]; kwargs...)
    end
    return Y
end

thetafilter(args...; pass = [2, 8], kwargs...) = bandpass(args...; pass, kwargs...)
gammafilter(args...; pass = [30, 400], kwargs...) = bandpass(args...; pass, kwargs...)
broadbandfilter(args...; pass = [10, 400], kwargs...) = bandpass(args...; pass, kwargs...)
narrowbandfilter(args...; pass = [50, 60], kwargs...) = bandpass(args...; pass, kwargs...)

# function alignlfp(x, y)
#     tx = dims(x, Ti)
#     y[Ti(Near(tx))]
# end

# harmonicbandpass(; kwargs...) = x -> harmonicbandpass(x; kwargs...)

# function harmonicbandpass(x::LFPVector; pass, harmonics=10)
#     t = dims(x, Ti)
#     T = t isa AbstractRange ? step(t) : t |> collect |> diff |> mean
#     fs = 1.0/T
#     _pass = [Tuple((i-1).*pass[1] .+ pass) for i in 1:harmonics]
#     display(_pass)
#     # ..................then also need the stop bands..............
#     y = filtfilt(remez(35, [p => 1 for p in _pass]; Hz=fs), collect(x))
#     return DimArray(y, dims(x))
# end

# function bandpass(X::LFPMatrix, dim=Dim{:channel}; kwargs...)
#     Y = similar(X)
#     Threads.@threads for x in dims(X, dim)
#         Y[dim(At(x))] .= bandpass(X[dim(At(x))]; kwargs...)
#     end
#     return Y
# end

# harmonicthetafilter(args...; pass=[4, 6], kwargs...) = harmonicbandpass(args...; pass, kwargs...)

"""
Detect time series with strong theta events in the first hald. We will call these 'theta events'
"""
# function bandthetafeature(x::AbstractVector, fs; pass=4..6)#, harmonics=4)
#     if pass isa Interval
#         pass = [pass]
#     end
#     # if harmonics > 0 # also count the harmonics of pass. Useful for highly nonlinear signals
#     #     pass = [(minimum(pass[1])*i)..(maximum(pass[1])*(i+1)) for i = 1:harmonics]
#     # end
#     x = x[1:round(Int, length(x)/2)]
#     𝐟 = rfftfreq(length(x), fs)
#     x̂ = rfft(x|>collect)
#     P = abs.(x̂).^2
#     idxs = [𝐟 .∈ (p,) for p in pass]
#     idxs = reduce((x, y)->x .| y, idxs)
#     return sum(P[idxs])/sum(P) # Proportion of power in the theta band
# end

"""
Detect time series with strong theta events using the automutual information
"""
function thetafeature(x::AbstractVector{<:Number}, fs = nothing; τ = 50, durprop = 0.5)
    # ! fs not implemented. τ in time steps
    @assert τ isa Integer
    x = x[1:round(Int, length(x) * durprop)]
    return selfmutualinfo(x, [τ]) |> first
end

function thetafeature(x::LFPVector; kwargs...)
    fs = 1.0 / step(dims(x, Ti))
    thetafeature(x, fs; kwargs...)
end

"""
Calculate a feature profile for each channel in each region
"""
function stimuluspartition(session, probeids, structures, stim; inbrain = 200,
                           times = nothing, epoch = nothing, kwargs...)
    Y = Vector{AbstractVector}([])
    stim == "spontaneous" &&
        return spontaneouspartition(session, probeids, structures; inbrain, kwargs...)
    epochs = selectepochs(session, stim, epoch)

    for p in eachindex(probeids)
        if isnothing(times)
            Z = []
            for e in axes(epochs, 1)
                _times = epochs[e, :].start_time .. epochs[e, :].stop_time
                if isnothing(structures)
                    X = getlfp(session, probeids[p]; inbrain, times = _times) |>
                        rectifytime
                else
                    X = getlfp(session, probeids[p], structures[p]; inbrain,
                               times = _times) |>
                        rectifytime
                end
                X = alignlfp(session, X, Symbol(epochs[e, :].stimulus_name); kwargs...)
                push!(Z, X)
            end
            X = vcat(Z...)
        else
            if isnothing(structures)
                X = getlfp(session, probeids[p]; inbrain, times) |>
                    rectifytime
            else
                X = getlfp(session, probeids[p], structures[p]; inbrain, times) |>
                    rectifytime
            end
            X = alignlfp(session, X, stim; kwargs...)
        end
        push!(Y, X)
    end
    return Y
end

function stimuluspartition(session, structures, stim; kwargs...)
    if isnothing(structures)
        probeids = getprobes(session).id
    else
        probeids = getprobes(session, structures)
    end
    stimuluspartition(session, probeids, structures, stim; kwargs...)
end

function spontaneouspartition(session, probeids, structures; inbrain = 200, mindur = 30,
                              kwargs...)
    times = stimulusintervals(session, "spontaneous").interval
    idxs = (diff.(extrema.(times) .|> collect) .|> first) .> mindur
    Y = Vector{AbstractVector}([])
    for p in eachindex(probeids)
        @info "Loading LFP for probe $p"
        X = []
        for (i, t) in enumerate(times)
            try
                if isnothing(structures)
                    _X = getlfp(session, probeids[p]; inbrain, times = t) |>
                         rectifytime
                else
                    _X = getlfp(session, probeids[p], structures[p]; inbrain, times = t) |>
                         rectifytime
                end
                push!(X, _X)
            catch e # Invalid interval, somehow
                @warn "Bad LFP for this probe at this time"
                _X = []
                push!(X, _X)
            end
        end
        push!(Y, X)
    end
    return Y
end

function spontaneouspartition(session, probeids, structures, duration; inbrain = 200)
    epoch = getepochs(session, "spontaneous")[2, :]
    times = epoch.start_time .. epoch.stop_time
    Y = Vector{AbstractVector}([])
    for p in eachindex(probeids)
        X = getlfp(session, probeids[p], structures[p]; inbrain, times) |> rectifytime
        _X = []
        t = minimum(dims(X, Ti))
        while t + duration < maximum(dims(X, Ti))
            push!(_X, X[Ti(t .. (t + duration))])
            t += duration
        end
        push!(Y, _X)
    end
    return Y
end

"""
Calculate the thetafeature for each stimulus presentation
"""
function thetafeature(Y::Vector{<:AbstractVector}; kwargs...) # Input formatted as stimuluspartition(
    t = [[([mean(dims(c, Ti)) for c in eachcol(s)]) for s in y] for y in Y]
    F = [[([thetafeature(c; kwargs...) for c in eachcol(s)]) for s in y] for y in Y]
    return F, t
end
function thetafeature(Y::Vector{<:AbstractVector{<:AbstractMatrix}}; kwargs...) # Input formatted as stimuluspartition(
    t = [[([mean(dims(c, Ti)) for c in eachcol(s)]) for s in y] for y in Y]
    F = [[([thetafeature(c; kwargs...) for c in eachcol(s)]) for s in y] for y in Y]
    return F, t
end
function thetafeature(Y::Vector{<:AbstractMatrix{<:Number}}; kwargs...) # Input formatted as stimuluspartition(
    t = [([mean(dims(c, Ti)) for c in eachcol(s)]) for s in Y]
    F = [([thetafeature(c; kwargs...) for c in eachcol(s)]) for s in Y]
    return F, t
end

# function meanfeature(Y::AbstractVector{AbstractVector{AbstractVector}})
#     [[([mean(dims(c, Ti)) for c in eachcol(s)]) for s in y] for y in Y]
# end

function ica(X::LFPMatrix, k = ceil(Int, size(X, 2) / 2))
    I = fit(ICA, collect(X)', k; maxiter = 1000)
    _X = DimArray(predict(I, collect(X)')', (dims(X, Ti), Dim{:channel}(1:size(I, 2))))
end

function pca(X::LFPMatrix)
    P = fit(PCA, collect(X)')
    _X = DimArray(predict(P, collect(X)')', (dims(X, Ti), Dim{:channel}(1:size(P, 2))))
    v = principalvars(P)
    return _X, v
end

DSP.hilbert(X::LFPMatrix) = mapslices(hilbert, X, dims = Ti)
DSP.hilbert(X::LFPVector) = DimArray(hilbert(X |> Array), dims(X); refdims = refdims(X))

function _waveletfreqs(t; moth = Morlet(2π), β = 1, Q = 32)
    n = length(t)
    fs = 1.0 ./ step(t) # Assume rectified time dim
    W = ContinuousWavelets.computeWavelets(n, wavelet(moth; β, Q);)[1]
    freqs = getMeanFreq(W, fs)
    freqs[1] = 0
    return freqs
end
function waveletfreqs(t; pass = nothing, kwargs...)
    freqs = _waveletfreqs(t; kwargs...)
    isnothing(pass) && return freqs
    return freqs[freqs .∈ [ClosedInterval(0, maximum(pass))]]
end

function _wavelettransform(x::AbstractVector; moth, β, Q) # β = 1 means linear in log space
    c = wavelet(moth; β, Q)
    res = ContinuousWavelets.cwt(x, c)
end

function _wavelettransform(t, x::AbstractVector; pass = nothing, moth, β, Q)
    if isnothing(pass)
        return _wavelettransform(x; moth, β, Q)
    end
    n = size(x, 1)
    @assert length(t) == n
    c = wavelet(moth; β, Q)
    W = ContinuousWavelets.computeWavelets(n, c)[1]
    freqs = getMeanFreq(W, 1.0 ./ step(t))
    pass = ClosedInterval(0, maximum(pass))
    W = W[:, freqs .∈ [pass]]
    res = ContinuousWavelets.cwt(x, c, W)[:, freqs .∈ [pass]]
end

function _wavelettransform(x::LFPVector; rectify = true, moth = Morlet(2π), β = 1, Q = 32,
                           pass = nothing)
    rectify && (x = rectifytime(x))
    t = dims(x, Ti)
    res = _wavelettransform(t, x |> Array; moth, β, Q, pass)
    freqs = waveletfreqs(t; moth, β, Q, pass)
    res = DimArray(res, (t, Dim{:frequency}(freqs)); metadata = DimensionalData.metadata(x),
                   refdims = refdims(x))
end

function _wavelettransform(x::LFPVector, ::Val{:mmap}; window = 50000, kwargs...)
    md = DimensionalData.metadata(x)
    rd = DimensionalData.refdims(x)
    x = rectifytime(x)
    𝓍 = _slidingwindow(x, window; tail = :overlap)
    t = dims(x, Ti)
    e = step(t) / 2
    freqs = waveletfreqs(dims(𝓍[1], Ti); kwargs...)
    sz = (length(t), length(freqs))
    fname = tempname()
    s = open(fname, "w+")
    write.((s,), sz)
    W = mmap(s, Matrix{ComplexF32}, sz)
    res = DimArray(W, (t, Dim{:frequency}(freqs)); metadata = (; md..., file = fname),
                   refdims = rd)
    threadlog, threadmax = (0, length(𝓍))
    @withprogress name="Wavelet transform" begin
        for _x in 𝓍
            subres = _wavelettransform(_x; rectify = false, kwargs...)
            tx = extrema(dims(subres, 1))
            fx = extrema(dims(subres, 2))
            tilims = Interval{:closed, :closed}(tx[1] - e, tx[2] + e)
            flims = Interval{:closed, :closed}(fx[1] - e, fx[2] + e)
            res[Ti(tilims), Dim{:frequency}(flims)] .= subres
            if threadmax > 1
                Threads.threadid() == 1 && (threadlog += 1) % 1 == 0 &&
                    @logprogress threadlog / threadmax
            end
        end
    end
    close(s)
    return res
end

_wavelettransform(x, s::Symbol; kwargs...) = _wavelettransform(x, Val(s); kwargs...)

function wavelettransform(x::LFPVector, args...; kwargs...)
    abs.(_wavelettransform(x, args...; kwargs...))
end

function fooof(p::LogWaveletMatrix, freqrange = [1.0, 300.0])
    ffreqs = 10.0 .^ collect(dims(p, Dim{:logfrequency}))
    freqrange = pylist([freqrange[1], freqrange[2]])
    spectrum = vec(collect(p))
    fm = PyFOOOF.FOOOF(peak_width_limits = pylist([0.5, 50.0]), max_n_peaks = 10,
                       aperiodic_mode = "knee", peak_threshold = 0.5)
    # if fm.aperiodic_params_[2] < 0.0
    #     fm = PyFOOOF.FOOOF(peak_width_limits=[0.5, 20.0], max_n_peaks=4, aperiodic_mode="fixed")
    # end
    # fm.report(freqs, spectrum, freqrange)
    fm.add_data(Py(ffreqs).to_numpy(), Py(spectrum).to_numpy(), freqrange)
    fm.fit()
    return fm
end

function logaperiodicfit(p::LogWaveletMatrix, freqrange = [1.0, 300.0], args...;
                         doplot = false)
    fm = fooof(p, freqrange, args...)
    if doplot != false
        p = fm.plot(; plt_log = true, file_name = doplot, save_fig = true)
    end
    # * The aperiodic model, as described in doi.org/10.1038/s41593-020-00744-x
    b, k, χ = pyconvert.((Float64,), fm.aperiodic_params_)
    k = max(k, 0.01)
    # b, k, χ = length(ps) == 3 ? ps : (ps[1], 0.0, ps[2])
    L = f -> 10.0 .^ (b - log10(k + (10.0^(f))^χ)) # Expects log frequency values
end

aperiodicfit(args...) = f -> (logaperiodicfit(args...)(log10(f)))

function _fooofedwavelet(res::LogWaveletMatrix, args...; kwargs...)
    psd = mean(res, dims = Ti)
    ffreqs = dims(psd, Dim{:logfrequency}) |> collect
    L = logaperiodicfit(psd, args...; kwargs...)
end
_fooofedwavelet(res::WaveletMatrix) = _fooofedwavelet(convert(LogWaveletMatrix, res))

function fooofedwavelet!(res::LogWaveletMatrix, freqrange = [1.0, 300.0]; kwargs...)
    psd = mean(res, dims = Ti)
    L = logaperiodicfit(psd, freqrange; kwargs...)
    ffreqs = dims(res, Dim{:logfrequency}) |> collect
    # f = Makie.Figure()
    # ax = Axis(f[1, 1]; xlabel="Frequency (Hz)")
    # Makie.lines!(ax, 10.0.^ffreqs[10:end], psd[:][10:end], color=:cornflowerblue)
    # Makie.lines!(ax, 10.0.^ffreqs[10:end], (L.(ffreqs)[10:end]), color=:crimson)
    # f
    # Makie.lines(ffreqs[10:end], psd[:][10:end].-L.(ffreqs)[10:end], color=:cornflowerblue)
    l = DimArray(L.(ffreqs), (Dim{:logfrequency}(ffreqs),))
    fmin = log10(minimum(freqrange))
    l[ffreqs .< fmin] = psd[ffreqs .< fmin] # Ensures no funky behavior outside of the fit bounds
    for r in axes(res, Ti)
        res[Ti(r)] .= res[Ti(r)] - l
    end
    return l
end

function fooofedwavelet(res::LogWaveletMatrix)
    _res = deepcopy(res)
    fooofedwavelet!(_res)
    return _res
end

fooofedwavelet(res::WaveletMatrix) = fooofedwavelet(convert(LogWaveletMatrix, res))

function fooofedwavelet(x::LFPVector; kwargs...)
    res = wavelettransform(x; kwargs...)
    fooofedwavelet(res)
end
import AllenNeuropixelsBase.PSDVector, AllenNeuropixelsBase.PSDMatrix
function aperiodicfit(psd::PSDVector, freqrange = [1.0, 300.0])
    ffreqs = dims(psd, Dim{:frequency}) |> collect
    freqrange = pylist([(freqrange[1]), (freqrange[2])])
    spectrum = vec(collect(psd))
    fm = PyFOOOF.FOOOF(peak_width_limits = pylist([0.5, 50.0]), max_n_peaks = 10,
                       aperiodic_mode = "knee", peak_threshold = 0.5)
    fm.add_data(Py(ffreqs).to_numpy(), Py(spectrum).to_numpy(), freqrange)
    fm.fit()
    b, k, χ = pyconvert.((Float64,), fm.aperiodic_params_)
    k = max(k, 0.01)
    L = f -> 10.0 .^ (b - log10(k + (f)^χ))
end

function fooofedspectrum!(psd::PSDMatrix, freqrange = [1.0, 300.0]; kwargs...)
    L = aperiodicfit.(eachcol(psd), (freqrange,); kwargs...)
    ffreqs = dims(psd, Dim{:frequency}) |> collect
    L = hcat([l.(ffreqs) for l in L]...)
    psd .= psd .- L
end

function wavelettransform(X::LFPMatrix; kwargs...)
    res = []
    threadlog, threadmax = (0, size(X, 2) / Threads.nthreads())
    @withprogress name="Wavelet transform" begin
        Threads.@threads for c in 1:axes(X, 2)
            push!(res, wavelettransform(X[:, c]; kwargs...))
            Threads.threadid() == 1 && (threadlog += 1) % 1 == 0 &&
                @logprogress threadlog / threadmax
        end
    end
    ti = dims(res[1], 1)
    freq = dims(res[1], 2)
    channel = dims(X, 2)
    return DimArray(cat(collect.(res)..., dims = 3), (ti, freq, channel))
end

mmapwavelettransform(x::LFPVector; kwargs...) = wavelettransform(x, :mmap; kwargs...)

function wavelettransform!(res::Dict, LFP::LFPMatrix; window = false, kwargs...)
    threadlog, threadmax = (0, size(LFP, 2) / Threads.nthreads())
    @withprogress name="Wavelet transform" begin
        for c in axes(LFP, 2)
            x = LFP[:, c]
            if window > 0
                # * Window the time series and calculate the wavelet transform in manageable chunks
                𝓍 = _slidingwindow(x, window)
                t = dims(x, Ti)
                freqs = waveletfreqs(dims(𝓍[1], Ti); kwargs...)
                sz = (length(t), length(freqs))
                s = open(tempname(), "w+")
                write.((s,), sz)
                W = mmap(s, Matrix{Float32}, sz)
                _res = DimArray(W, (t, Dim{:frequency}(freqs)))
                for _x in 𝓍
                    subres = wavelettransform(_x; kwargs...)
                    tx = extrema(dims(subres, 1))
                    fx = extrema(dims(subres, 2))
                    tilims = Interval{:closed, :closed}(tx[1] - eps(), tx[2] + eps())
                    flims = Interval{:closed, :closed}(fx[1] - eps(), fx[2] + eps())
                    _res[Ti(tilims), Dim{:frequency}(flims)] .= subres
                end
            else
                push!(res, dims(LFP, 2)[c] => wavelettransform(x; kwargs...)) # Doesnt actually write to file
            end
            Threads.threadid() == 1 && (threadlog += 1) % 1 == 0 &&
                @logprogress threadlog / threadmax
        end
    end
end

function TimeseriesSurrogates.surrogate(x::LFPVector, S::Surrogate; kwargs...)
    (y = deepcopy(x); y .= surrogate(x |> collect .|> Float64, S; kwargs...) .|> eltype(y); y)
end

function TimeseriesSurrogates.surrogenerator(x::LFPVector, S::IAAFT)
    sg = surrogenerator(x |> collect .|> Float64, S)
    return () -> DimArray(sg() .|> eltype(x), dims(x))
end
function TimeseriesSurrogates.surrogenerator(X::LFPMatrix, S::IAAFT)
    sg = [surrogenerator(x |> collect .|> Float64, S) for x in eachcol(X)]
    function out()
        Y = deepcopy(X)
        [y .= s() .|> eltype(X) for (s, y) in zip(sg, eachcol(Y))]
        return Y
    end
    return out
end

function powerspectra(t, X; n = 5000, window = DSP.Windows.hanning)
    Δt = t[2] - t[1]
    @assert all(Δt .≈ diff(t))
    fp = x -> welch_pgram(x, n; fs = 1 / Δt, window)
    P = [fp(Array(x)) for x in eachcol(X)]
    𝑓 = P[1].freq
    psd = hcat([p.power for p in P]...)
    # psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
    return 𝑓, psd
end

function powerspectra(t, x, X; kwargs...)
    𝑓, psd = powerspectra(t, X; kwargs...)
    psd = DimArray(psd, (Dim{:frequency}(𝑓), Dim{:channel}(x)))
end

function powerspectra(X::LFPMatrix; kwargs...)
    t = dims(X, Ti) |> collect
    x = dims(X, Dim{:channel}) |> collect
    return powerspectra(t, x, X; kwargs...)
end

function powerspectra(x::LFPVector; kwargs...)
    t = dims(x, Ti) |> collect
    𝑓, psd = powerspectra(t, collect(x); kwargs...)
    psd = DimArray(psd[:], (Dim{:frequency}(𝑓),))
end

function mutualinfo(x, y, est; base = 2, α = 1)
    X = genentropy(Dataset(x), est; base, α)
    Y = genentropy(Dataset(y), est; base, α)
    XY = genentropy(Dataset(x, y), est; base, α)
    return X + Y - XY
end

function phasematch(a::Tuple, b::Tuple; tol = 0.05) # a = (LFP1, PHI1)
    x, xp = a
    y, yp = b
    i = findfirst(0 .< (yp .- xp[end]) .< tol) # Phase is always increasing
    if isnothing(i)
        _, ix = findmin(abs.(yp .- xp[end]))
        ss = (yp .- xp[end])
        @warn "No matching phases found. b is of length $(length(b[1])). The final phase of a is $(xp[end]). The extrema of phases in b is $(extrema(yp)). The smallest difference in phase is $(ss[ix]), at index $ix."
        return a
    else
        x = vcat(x, y[i:end])
        xp = vcat(xp, yp[i:end])
    end
    return (x, xp)
end

function phasematch(a::AbstractVector, b::AbstractVector, pass = nothing; kwargs...)
    if !isnothing(pass)
        _a = bandpass(deepcopy(a); pass)
        _b = bandpass(deepcopy(b); pass)
    else
        _a = deepcopy(a)
        _b = deepcopy(b)
    end
    pha = hilbert(_a) .|> angle
    phb = hilbert(_b) .|> angle
    return phasematch((a, pha), (b, phb); kwargs...)
end

function phasematch(a::Vector, b::Vector, pass = nothing; kwargs...)
    if !isnothing(pass)
        _a = bandpass(deepcopy(a); pass)
        _b = bandpass(deepcopy(b); pass)
    else
        _a = deepcopy(a)
        _b = deepcopy(b)
    end
    pha = hilbert(_a) .|> angle
    phb = hilbert(_b) .|> angle
    return phasematch((a, pha), (b, phb); kwargs...)
end

# phasematch(ab::Tuple; kwargs...) = phasematch(ab[1], ab[2], pass=nothing; kwargs...)

function _extracttheta(session, stimulus, structures; inbrain = 200, times = nothing,
                       trail = :offset, thetathresh = 0.6, durprop = 0.5)
    probeids = getprobes(session, structures)

    Y = stimuluspartition(session, probeids, structures, stimulus; inbrain, times, trail)

    τ = 50 # Approx min of global AMI

    durprop # Events tend to last for half of the flashes stimulus interval. Set this to only take the first 0.5 of the presentation TS

    F, t = thetafeature(Y; τ, durprop)
    F̄ = [[mean(x) for x in X] for X in F]
    t̄ = [[mean(x) for x in X] for X in t]

    # Extract the presentations with a thetafeature above 0.5 in VISp
    idxs = F̄[1] .> thetathresh
    f = [F̄[i][idxs] for i in eachindex(F̄)]
    Y = [Y[i][idxs] for i in eachindex(Y)]
    Y = [[x[1:round(Int, size(x, 1) * durprop), :] for x in y] for y in Y]
    t = [t̄[i][idxs] for i in eachindex(t̄)]

    return Y, f
end
function _extracttheta(session::Int, args...; kwargs...)
    _extracttheta(Session(session), args...; kwargs...)
end

function extracttheta(session, stimulus, structures; cattimes = false, kwargs...)
    Y, f = _extracttheta(session, stimulus, structures; kwargs...)
    # Y = reduce.((phasematch,), Y; pass=[1, 10])
    if cattimes
        Y = [cat(y...; dims = Ti) for y in Y]
    else
        Y = catlfp.(Y)
    end
end
function extracttheta(session::Int, args...; kwargs...)
    extracttheta(Session(session), args...; kwargs...)
end
function extracttheta(params::NamedTuple, args...; kwargs...) # Assume the pair is "VISp", "VISl"
    if params[:stimulus] == "flashes"
        structures = ["VISp", "VISl"]
        x = extracttheta(params[:sessionid], params[:stimulus], structures, args...;
                         kwargs...)
        structure = params[:structure]
        idx = findfirst(structures .== structure)
        return x[idx]
    else
        return formatlfp(; params...)
    end
end
function _extracttheta(params::NamedTuple, args...; kwargs...) # Assume the pair is "VISp", "VISl"
    if params[:stimulus] == "flashes"
        structures = ["VISp", "VISl"]
        x, _ = _extracttheta(params[:sessionid], params[:stimulus], structures, args...;
                             kwargs...)
        structure = params[:structure]
        idx = findfirst(structures .== structure)
        return x[idx]
    else
        return formatlfp(; params...)
    end
end
