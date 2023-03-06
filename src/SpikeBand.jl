using SparseArrays
using ProgressLogging
using Distances
using Clustering
using Random

function downloadspikes(S::AbstractSession)
    _ = S.pyObject.spike_times
    _ = S.pyObject.spike_amplitudes
    return nothing
end

SpikeMatrix = SparseDimArray{T, 2, Tuple{A, B}} where {T, A<:DimensionalData.TimeDim, B<:Dim{:unit}}
export SpikeMatrix

function getsessionpath(session::AbstractSession)
    path = joinpath(datadir, "Ecephys", "session_"*string(getid(session)), "session_"*string(getid(session))*".nwb");
end

getspiketimes(S::AbstractSession) = S.pyObject.spike_times
getspikeamplitudes(S::AbstractSession) = S.pyObject.spike_amplitudes

function getspiketimes(S::AbstractSession, structure::AbstractString)
    unitstructs = getunitmetrics(S)
    unitids = unitstructs[unitstructs.structure_acronym .== structure, :].unit_id
    spiketimes = filter(p -> p[1] in unitids, getspiketimes(S))
end

function getspiketimes(S::AbstractSession, unitids::Vector{Int})
    spiketimes = getspiketimes(S)
    Dict(a => b for (a,b) in spiketimes if a in unitids)
end

function formatspiketimes(; sessionid=757216464, structure="VISp", stimulus="gabors", n=1, filter=true, kwargs...)
    S = Session(sessionid)
    spiketimes = getspiketimes(S, structure)
    I = stimulusepochs(S, stimulus).interval[n]
    ks = keys(spiketimes)
    vals = values(spiketimes)
    vals = [v[in.(v, (I,))] for v in vals]
    if filter
        metrics = getunitanalysismetrics(S)
        Sp = Dict(k => v for (k,v) in zip(ks, vals) if k ∈ metrics.ecephys_unit_id)
    else
        Sp = Dict(k => v for (k,v) in zip(ks, vals))
    end
end


"""
Construct a sparse array of spike counts from a Dict of spike times
"""
function spikematrix(Sp::AbstractDict, bin=1e-4; rectify=4)
    tmin, tmax = extrema(vcat(values(Sp)...))
    units = Sp |> keys |> collect
    if rectify > 0
        tmax = round(tmax; digits=rectify)
        tmin = round(tmin; digits=rectify)
    end
    ts = tmin:bin:tmax
    spikes = spzeros(Float32, length(ts), length(units))
    spikes = goSparseDimArray(spikes, (Ti(ts), Dim{:unit}(units)))
    @withprogress name="spikearray" begin
        for u in eachindex(units)
            _times = Sp[units[u]]
            for t in eachindex(_times) # Hella slow
                spikes[Ti(Near(_times[t])), Dim{:unit}(At(units[u]))] += 1
            end
            @logprogress u/length(units)
        end
    end
    return spikes
end

function burstmatrix(Bs::AbstractDict, bin=1e-2; rectify=4)
    cs = collect(keys(Bs))
    b = spikematrix(Bs, bin; rectify)
    b = goSparseDimArray(b.data, (Ti(collect(dims(b, Ti))), Dim{:channel}(cs)))
end
function burstmatrix(B::AbstractVector{<:BurstVector}, args...; kwargs...)
    channels = [getchannel(b[1]) for b in B]
    d = Dict(channels[i] => peaktime.(B[i]) for i in eachindex(channels))
    return burstmatrix(d, args...; kwargs...)
end

function smoothspikematrix(s; windowfunc=hanning, window=0.1)
    ts = collect(s.dims[1])
    dt = mean(diff(collect(s.dims[1])))
    n = ceil(Int, window/dt/2)*2
    x = hcat([DSP.conv(x[:], windowfunc(n)) for x in eachcol(collect(s))]...)
    x = x[n÷2:end-n÷2, :]
    return rectifytime(DimArray(x, (Ti(ts), dims(s, 2))))
    # ! Note that the hanning window has an integral of a half

end


# `count` is a boolean indicating whether to return the number of spikes in each bin, or sum the amplitudes
function _getspikes(units, times, amplitudes, _times, bin, rectify, count)
    tmax = maximum(_times)
    tmin = minimum(_times)
    times = [filter(∈(tmin..tmax), t) for t in times]
    if rectify > 0
        tmax = round(tmax; digits=rectify)
        tmin = round(tmin; digits=rectify)
    end
    ts = tmin:bin:tmax
    spikes = spzeros(Float32, length(ts), length(units))
    spikes = goSparseDimArray(spikes, (Ti(ts), Dim{:unit}(units))) # SparseDimArray?
    @withprogress name="spikearray" begin
    for u in eachindex(units)
        _times = times[u]
        _amplitudes = amplitudes[u]
        for t in eachindex(_times) # Hella slow
            spikes[Ti(Near(_times[t])), Dim{:unit}(At(units[u]))] += (count ? 1 : _amplitudes[t])
        end
        @logprogress u/length(units)
    end
    end
    return spikes
end

"""
We combine the spike times and spike amplitudes into one sparse array, using a given bin width.
"""
function getspikes(S, timebounds=nothing; bin=1e-4, rectify=4, structure=nothing, count=true)
    times = isnothing(structure) ? getspiketimes(S) : getspiketimes(S, structure)
    validunits = findvalidunits(S, keys(times))
    times = Dict(k => v for (k,v) in times if k ∈ validunits)
    units = times |> keys |> collect
    times = times |> values |> collect
    amplitudes = S |> getspikeamplitudes
    amplitudes = getindex.((amplitudes,), units)
    @assert all(length.(times) .== length.(amplitudes))
    # Set up the sparse dim array, then allocate spikes to the nearest time index
    _times = isnothing(times) ? vcat(_times...) : timebounds
    _getspikes(units, times, amplitudes, _times, bin, rectify, count)
end

function getspikes(S, stimulus::String; n=1, kwargs...)
    timebounds = getstimulustimes(S, stimulus)[n]
    getspikes(S, timebounds; kwargs...)
end

function getspikes(S, stimulus::String, structure::String; kwargs...)
    spikes = getspikes(S, stimulus; kwargs...)
    structures = getstructureacronyms(S, dims(spikes, Dim{:unit}))
    spikes = spikes[:, findall(structures .== structure)]
end

"""
A function to easily grab formatted spike data for a given session, using some sensible default parameters
"""
function formatspikes(; sessionid=757216464, stimulus="gabors", structure="VISp", epoch = 1, bin=1e-4, kwargs...)
    sesh = Session(sessionid)
    S = getspikes(sesh, stimulus, structure; bin, rectify=4, structure, count=true, epoch)
end



function getstructureacronyms(session::AbstractSession, units)
    unittable = getunitmetrics(session)
    acronyms = Vector{Any}(undef, size(units))
    [acronyms[i] = notemptyfirst(unittable[unittable.unit_id.==units[i], :ecephys_structure_acronym]) for i ∈ 1:length(units)]
    return acronyms
end


function minspikediff(S::AbstractSession)
    diffs = S |> getspiketimes |> values .|> diff
    return minimum(minimum.(diffs)) # Just greater than 1e-4
end


# function clusterunits(spikes; dist=CosineDist())
#     D = pairwise(dist, spikes|>Array, dims=2)
#     h = hclust(D)
# end


function getreceptivefield(session, unit)
    rf = stimulusmapping.ReceptiveFieldMapping(session.pyObject)
    field = rf.get_receptive_field(unit)
end

function getisis(x::AbstractSparseDimVector)
    s = findnz(x |> SparseVector)[1]
    t = dims(x, Ti) |> collect
    I = AN.ephys_features.get_isis(t, s)
    @assert I ≈ diff(t[s]) # Yeah there's this
    return I
end

function getisis(S::AbstractSession)
    ts = getspiketimes(S)
    isis = Dict(a => diff(b) for (a,b) in ts)
end

function getisis(S::AbstractSession, units)
    isis = getisis(S)
    isis = Dict(a => diff(b) for (a,b) in isis if a in units)
end

# function detectbursts(isis::AbstractVector)
#     AN.ephys_features.detect_bursts(isis, isi_types, fast_tr_v, fast_tr_t, slow_tr_v, slow_tr_t, thr_v, tol=0.5, pause_cost=1.0)
# end

import DataFrames.subset
function subset(d::DataFrame, col, vals::AbstractVector)
    idxs = indexin(vals, d[:, col])
    return d[idxs, :]
end
subset(d::DataFrame, col, vals::Base.KeySet) = subset(d, col, collect(vals))
subset(d::DataFrame, col, vals::Dim) = subset(d, col, collect(vals))
subset(d::DataFrame, col, vals::DataFrame) = subset(d, col, vals[:, col])
function subset(d::DataFrame, col, val)
    idxs = d[:, col] .== val
    return d[idxs, :]
end
export subset


# function receptivefieldcentrality(unit; centre=(0, 0), metrics=AN.getunitanalysismetricsbysessiontype("brain_observatory_1.1"))
#     # rf = AN.getreceptivefield(session, unit)
# end

function receptivefieldfilter(am::DataFrame)
    am = am[[ismissing(𝑝) || 𝑝 < 0.01 for 𝑝 in am.p_value_rf], :]
    return am.unit_ids
end

function receptivefieldfilter(units::AbstractVector; am = AN.getunitanalysismetricsbysessiontype("brain_observatory_1.1"))
    am = AN.subset(am, :ecephys_unit_id, units)
    return receptivefieldfilter(am)
end



function alignspiketimes(session, X, ::Val{:flashes}; x_position=nothing, y_position=nothing)
    stims = stimulusintervals(session, "flashes")
    intervals = stims.interval
    times = values(X)
    units = keys(X)
    _times = []
    for ts in times
        _ts = []
        for i in intervals
            push!(_ts, ts[ts .∈ (i,)])
        end
        push!(_times, _ts)
    end
    return Dict(units .=> _times)
end

alignspiketimes(session, X, stimulus="flashes"; kwargs...) = alignspiketimes(session, X, stimulus|>Symbol|>Val; kwargs...)





function countspikes(ts::AbstractVector, T::Real)
    # ts should be sorted
    ts = deepcopy(ts)
    m = maximum(ts)
    _t = minimum(ts)
    x = zeros(floor(Int, (m-_t)/T))
    for i in eachindex(x)
        idxs = (ts .≥ (_t + (i-1)*T)) .& (ts .< (_t + i*T,))
        x[i] = length(splice!(ts, findall(idxs))) # * Make x shorter so future searches are quicker
    end
    return x
end

function fanofactor(ts::AbstractVector, T::Real)
    # * Calculate the fano factor using non-overlapping windows
    ts = sort(ts)
    x = countspikes(ts, T)
    return var(x)/mean(x)
end

function defaultfanobins(ts)
    maxwidth = (first∘diff∘collect∘extrema)(ts)/10
    minwidth = max((mean∘diff)(ts), maxwidth/10000)
    # spacing = minwidth
    return 10.0.^range(log10(minwidth), log10(maxwidth); length=100)
end

fanofactor(ts::AbstractVector, T::AbstractVector=defaultfanobins(ts)) = (T, fanofactor.((ts,), T))




function metricmap(stimulus)
    if stimulus == "flashes"
        return :fl
    elseif stimulus == "drifting_gratings"
        return :dg
    elseif stimulus == "natural_scenes"
        return :ns
    elseif stimulus == "gabors"
        return :rf
    end
end

function getmetric(metrics::DataFrame, metric, stimulus)
    sesh = metricmap(stimulus)
    return metrics[:, Symbol(reduce(*, string.([metric, :_, sesh])))]
end

function getfano(metrics::DataFrame, stimulus)
    sesh = metricmap(stimulus)
    return metrics[:, Symbol(reduce(*, string.([:fano_, sesh])))]
end

function findvalidunits(session, units; kwargs...)
    units = collect(units)
    metrics = getunitanalysismetrics(session; filter_by_validity=true, kwargs...)
    check = units .∈ (metrics.ecephys_unit_id,)
    return units[check]
end

function findvalidunits(session, spikes::Dict)
    units = findvalidunits(session, keys(spikes))
    return Dict(units .=> getindex.((spikes,), units))
end


function lfpspikesimilarity(X::LFPVector, Y::Vector; normalize=identity)
    ts = Interval(extrema(dims(X, Ti))...)
    Y = Y[Y .∈ (ts,)]
    X = normalize(X)
    return mean(X[Ti(Near(Y))]) # Mean value of the LFP at each spike time
end

function lfpspikesimilarity(session, probeid, X::LFPMatrix, Y::Dict; normalize=identity, kwargs...)
    X = normalize(X)
    spikes = values(Y) |> collect
    units = keys(Y) |> collect
    # * First, match up the units to their nearest LFPs
    depths = getchanneldepths(session, probeid, X)
    idxs = sortperm(depths)
    X = X[:, idxs]
    depths = depths[idxs]
    X = DimArray(X, (dims(X, Ti), Dim{:depth}(depths)))
    unitdepths = getunitdepths(session, probeid, units)
    X = X[Dim{:depth}(Near(unitdepths))] # Sorted in order of spikes
    sims = [lfpspikesimilarity(X[:, i], spikes[i]; kwargs...) for i in eachindex(spikes)]
    return Dict(units .=> sims)
end

function getclosestchannels(session, probeid, units, channels)
    channels = collect(channels)
    cdepths = getchanneldepths(session, probeid, channels)
    udepths = getunitdepths(session, probeid, units)
    dists = [abs.(cdepths .- u) for u in udepths]
    idxs = last.(findmin.(dists))
    channels = channels[idxs]
    Dict(units .=> channels)
end
getclosestchannels(sessionid::Int, args...) = getclosestchannels(Session(sessionid), args...)

function getpeakchannels(session, units)
    metrics = getunitmetrics(session)
    check = all(units .∈ (metrics.unit_id,))
    if !check
        @error "Some units do not have corresponding metrics"
    end
    metrics = subset(metrics, :unit_id, units)
    channels = Dict(zip(units, metrics.peak_channel_id))
end
getpeakchannels(sessionid::Int, args...) = getpeakchannels(Session(sessionid), args...)


function getunitchannels(session, units)
    metrics = getunitanalysismetrics(session)
    check = all(units .∈ (metrics.ecephys_unit_id,))
    if !check
        @error "Some units do not have corresponding metrics"
    end
    metrics = subset(metrics, :ecephys_unit_id, units)
    channels = Dict(zip(metrics.ecephys_channel_id))
end
getunitchannels(sessionid::Int, args...) = getunitchannels(Session(sessionid), args...)




mape(x, y) = sum(1.0.-abs.((x .- y)./x))/length(x)
sse(x, y) = sum((x .- y).^2)
msse(x, y) = std((x .- y).^2)
R²(x, y) = 1.0 - sse(x, y)/sse(x, mean(x))
absoluteR²(x, y) = R²(x, y) > 0 ? R²(x, y) : R²(x, -y) # If the predictions are negatively correlated to x

function predictionerror(x, y, M::CCA; metric=R²)

    # * Forward prediction, x -> y
    zx = predict(M, x, :x)
    zy = zx # ? The assumption that the low-dim shared space maximises correlation between the two latent variables
    ŷ = (M.yproj)'\zy .+ M.ymean

    # * Reverse prediction, y -> x
    zy = predict(M, y, :y)
    zx = zy
    x̂ = (M.xproj)'\zx .+ M.xmean
    return metric(x, x̂), metric(y, ŷ)
end


"""
Calculate the prediction error of variables 1 to variables 2 and vice versa.
The output contains prediction errors as (x->y, y->x).
If used with spike matrices probably want to transpose those
"""
function predictionerror(x, y; metric=R², model=MultivariateStats.CCA, kwargs...)
    x = collect(x)
    y = collect(y)
    N = size(x, 2)
    nₜ = N÷5
    iₜ = fill(false, N)
    iₜ[randperm(N)[1:nₜ]] .= true
    xₜ = x[:, iₜ]
    yₜ = y[:, iₜ]
    x = x[:, .!iₜ]
    y = y[:, .!iₜ]
    N = min(size(x, 1), size(y, 1))
    Δ = [predictionerror(x, y, fit(model, x, y; outdim=n, kwargs...); metric) for n in 1:N]
    display(Δ)
    return first.(Δ), last.(Δ)
end
