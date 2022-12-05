using DimensionalData
using PyCall
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

LFPVector = AbstractDimArray{T, 1, Tuple{A}, B} where {T, A<:DimensionalData.TimeDim, B}
LFPMatrix = AbstractDimArray{T, 2, Tuple{A, B}} where {T, A<:DimensionalData.TimeDim, B<:Dim{:channel}}
export LFPMatrix, LFPVector # For simpler dispatch
dimmatrix(a::Symbol, b::Symbol) = AbstractDimArray{T, 2, Tuple{A, B}} where {T, A<:Dim{a}, B<:Dim{b}}
dimmatrix(a, b::Symbol) = AbstractDimArray{T, 2, Tuple{A, B}} where {T, A<:a, B<:Dim{b}}
dimmatrix(a, b) = AbstractDimArray{T, 2, Tuple{A, B}} where {T, A<:a, B<:b}
export dimmatrix

PSDMatrix = dimmatrix(:frequency, :channel)
PSDVector = AbstractDimArray{T, 1, Tuple{A}, B} where {T, A<:Dim{:frequency}, B}
LogPSDVector = AbstractDimArray{T, 1, Tuple{A}, B} where {T, A<:Dim{:logfrequency}, B}

duration(X::AbstractDimArray) = diff(extrema(dims(X, Ti))|>collect)|>first
function samplingperiod(X::AbstractDimArray)
    if dims(X, Ti).val.data isa AbstractRange
        step(dims(X, Ti))
    else
        mean(diff(collect(dims(X, Ti))))
    end
end
samplingrate(X::AbstractDimArray) = 1/samplingperiod(X)


WaveletMatrix = dimmatrix(Ti, :frequency) # Type for DimArrays containing wavelet transform info
LogWaveletMatrix = dimmatrix(Ti, :logfrequency) # Type for DimArrays containing wavelet transform info
export WaveletMatrix, LogWaveletMatrix
function Base.convert(::Type{LogWaveletMatrix}, x::WaveletMatrix)
    x = DimArray(x, (dims(x, Ti), Dim{:logfrequency}(log10.(dims(x, :frequency)))); metadata=metadata(x), refdims=refdims(x))
    x = x[:, .!isinf.(dims(x, :logfrequency))]
end
function Base.convert(::Type{LogPSDVector}, x::PSDVector)
    x = DimArray(x, (Dim{:logfrequency}(log10.(dims(x, :frequency))),); metadata=metadata(x), refdims=refdims(x))
    x = x[.!isinf.(dims(x, :logfrequency))]
end
function Base.convert(::Type{WaveletMatrix}, x::LogWaveletMatrix)
    x = DimArray(x, (dims(x, Ti), Dim{:frequency}(exp10.(dims(x, :logfrequency)))); metadata=metadata(x), refdims=refdims(x))
end
waveletmatrix(res::LogWaveletMatrix) = convert(WaveletMatrix, res)
logwaveletmatrix(res::WaveletMatrix) = convert(LogWaveletMatrix, res)
logpsdvector(res::PSDVector) = convert(LogPSDVector, res)


function downloadlfp(S::AbstractSession, probeid::Int)
    @assert(any(getprobeids(S) .== probeid), "Probe $probeid does not belong to session $(getid(S))")
    @assert(subset(getprobes(S), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    _ = S.pyObject.get_lfp(probeid)
    return nothing
end


function structure2probe(S::AbstractSession, structure::String)
    channels = getchannels(S)
    channels = channels[.!ismissing.(channels.structure_acronym), :]
    channels = subset(channels, :structure_acronym, structure)
    @assert length(unique(channels.probe_id)) == 1
    return channels.probe_id |> unique |> first
end

function getlfppath(session::AbstractSession, probeid)
    path = joinpath(datadir, "Ecephys", "session_"*string(getid(session)), "probe_"*string(probeid)*"_lfp.nwb");
end

function getlfptimes(session::AbstractSession, probeid)
    path = getlfppath(session, probeid)
    if !isfile(path)
        @error "Probe lfp data file does not exist. This may be becuase you have not downloaded the probe data. Try using `getlfp`"
    end

    f = h5open(path)
    times = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["timestamps"][:]
    close(f)
    return times
end

function getlfpchannels(session::AbstractSession, probeid)
    if !isfile(getlfppath(session, probeid))
        downloadlfp(session, probeid)
    end
    f = h5open(getlfppath(session, probeid))
    channels = f["general"]["extracellular_ephys"]["electrodes"]["id"][:]
    close(f)
    return channels
end


"""
Get the lfp data for a probe,  providing *indices* for channels and times. See function below for indexing by channel ids and time values/intervals
"""
function _getlfp(session::AbstractSession, probeid::Int; channelidxs=1:length(getlfpchannels(session, probeid)), timeidxs=1:length(getlfptimes(session, probeid)))
    @assert(any(getprobeids(session) .== probeid), "Probe $probeid does not belong to session $(getid(session))")
    @assert(subset(getprobes(session), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    path = getlfppath(session, probeid)
    if !isfile(path)
        downloadlfp(session, probeid)
    end
    if !isfile(path)
        @error "LFP data did not download correctly"
    end

    f = h5open(path)

    timedata = getlfptimes(session, probeid)
    timedata = timedata[timeidxs]
    dopermute = true
    channelids = getlfpchannels(session, probeid)
    channelids = channelids[channelidxs]
    if (channelidxs isa Union{Int64, AbstractRange{Int64}}) & (timeidxs isa Union{Int64, AbstractRange{Int64}}) # Can use HDF5 slicing
        lfp = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][channelidxs, timeidxs]
    elseif timeidxs isa Union{Int64, AbstractRange{Int64}}
        lfp = [f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][i, timeidxs] for i ∈ channelidxs]
        lfp = hcat(lfp...)
        dopermute = false
    else
        lfp = read(f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"])
        lfp = lfp[channelidxs, timeidxs]
    end
    if lfp isa Vector
       lfp = reshape(lfp, 1, length(lfp))
    end
    if channelids isa Number
        channelids = [channelids]
    end
    if dopermute
        lfp = permutedims(lfp, reverse(1:ndims(lfp)))
    end
    X = DimArray(lfp, (Ti(timedata),  Dim{:channel}(channelids)))
    close(f)
    return X
end


function isinvalidtime(session, probeids=getprobeids(session), times=NaN)
    if isempty(session.pyObject.get_invalid_times()) # No invalid times in this session!
        return false
    end
    intervals = [session.pyObject.get_invalid_times().start_time.values, session.pyObject.get_invalid_times().stop_time.values]
    intervals = [[intervals[1][i], intervals[2][i]] for i ∈ 1:length(intervals[1])]
    tags = session.pyObject.get_invalid_times().tags.values
    tags = vcat([tryparse.(Int64, i) for i ∈ tags]...)
    badprobes = tags[.!isnothing.(tags)]
    if times isa Interval
        isininterval = [any(i .∈ (times,)) for i ∈ intervals]
    else
        isininterval = [any((times .> i[1]) .& (times .< i[2])) for i ∈ intervals]
    end
    return any(probeids .∈ (badprobes,)) & any(isininterval)
end

"""
This is the one you should be using. Get lfp data by channel id and time intervals or vector. Also, throw error if you try to access an invalid time interval.
"""
function getlfp(session::AbstractSession, probeid::Int; channels=getlfpchannels(session, probeid), times=ClosedInterval(extrema(getlfptimes(session, probeid))...), inbrain=false)
    if isinvalidtime(session, probeid, times)
        @error "Requested LFP data contains an invalid time..."
    end
    if inbrain isa Symbol || inbrain isa Real || inbrain
        depths = getchanneldepths(session, probeid, channels)
        if inbrain isa Real # A depth cutoff
            channels = channels[depths .> inbrain]
        elseif inbrain isa Symbol # A mode
        else # Just cutoff at the surface
            channels = channels[depths .> 0]
        end
    end

    channelidxs = getlfpchannels(session, probeid)
    channelidxs = indexin(channels, channelidxs)
    channelidxs = filter(!isnothing, channelidxs)
    @assert length(channels) == length(channelidxs)

    timeidxs = getlfptimes(session, probeid)
    if !(times isa Interval) && length(times) == 2
        times = ClosedInterval(times...)
    end
    if times isa Interval
        timeidxs = findall(timeidxs .∈ (times,))
    else
        timeidxs = indexin(times, timeidxs)
        timeidxs = filter(!isnothing, timeidxs)
        @assert length(timeidxs) == length(times)
    end

    # See if we can convert to unitranges for faster HDF5 reading via slicing
    if collect(UnitRange(extrema(timeidxs)...)) == timeidxs
        timeidxs = UnitRange(extrema(timeidxs)...)
    end
    if collect(UnitRange(extrema(channelidxs)...)) == channelidxs
        channelidxs = UnitRange(extrema(channelidxs)...)
    end

    @info "Accessing LFP data"
    _getlfp(session, probeid; channelidxs, timeidxs)
end

"""
If you want to downsample the LFP data, its quicker to use this function and then perform slicing afterwards (since getlfp() has to check all of the time coordinates you supply, which can be slow).
"""
function getdownsampledlfp(session, probeid; downsample=100, timerange=ClosedInterval(extrema(getlfptimes(session, probeid))...), channels=getlfpchannels(session, probeid))
    if !(timerange isa Interval) && length(timerange) == 2
        timerange = ClosedInterval(timerange...)
    end
    timevals = getlfptimes(session, probeid)
    tidxs = timevals .∈ (timerange,)
    times = findfirst(tidxs):downsample:findlast(tidxs)
    _getlfp(session, probeid; timeidxs = times)[:, At(channels)]
end


"""
Now we can overload `getlfp()` to index by structure
"""
function getlfp(session::AbstractSession, probeid::Int, structures::Union{Vector{String}, String}; kwargs...)
    if structures isa String
        structures = [structures]
    end
    channels = subset(getchannels(session, probeid), :ecephys_structure_acronym=>ByRow(∈(structures)), skipmissing=true)
    channels = channels.id ∩ getlfpchannels(session, probeid)
    getlfp(session, probeid; channels, kwargs...)
end

function getlfp(S::AbstractSession, structure::String; kwargs...)
    probeid = structure2probe(S, structure)
    getlfp(S, probeid, structure; kwargs...)
end

function getlfp(session, probeids::Vector{Int}, args...; kwargs...)
    LFP = [getlfp(session, probeid, args...; kwargs...) for probeid ∈ probeids]
end

function formatlfp(; sessionid=757216464, probeid=769322749, stimulus="gabors", structure="VISp", n=1, kwargs...)
    session = Session(sessionid)
    if isnothing(structure)
        structure = getchannels(session, probeid).id |> getstructureacronyms |> unique |> skipmissing |> collect |> Vector{String}
    end
    epoch = getepochs(session, stimulus)[n, :]
    times = epoch.start_time..epoch.stop_time
    X = getlfp(session, probeid, structure; inbrain=200, times) |> rectifytime
    X = sortbydepth(session, probeid, X)
    X = DimArray(X; metadata=Dict(:sessionid=>sessionid, :probeid=>probeid, :stimulus=>stimulus, :structure=>structure))
end





function getchannels(data::AbstractDimArray)
    dims(data, :channel).val
end

"""
At the moment this is just a proxy: distance along the probe to the cortical surface
"""
function getchanneldepths(session, probeid, channels)
    cdf = getchannels(session, probeid)
    return _getchanneldepths(cdf, channels)
end
# function getchanneldepths(session, channels)
#     cdf = getchannels(session) # Slightly slower than the above
#     #cdf = cdf[indexin(channels, cdf.id), :]
#     cdfs = groupby(cdf, :probe_id)
#     depths = vcat([_getchanneldepths(c, c.id) for c ∈ cdfs]...)
#     depths = depths[indexin(channels, vcat(cdfs...).id)]
# end
function _getchanneldepths(cdf, channels)
    surfaceposition = minimum(subset(cdf, :ecephys_structure_acronym=>ByRow(ismissing)).probe_vertical_position)
    # Assume the first `missing` channel corresponds to the surfaceprobe_vertical_position
    idxs = indexin(channels, cdf.id)[:]
    alldepths = surfaceposition .- cdf.probe_vertical_position # in μm
    depths = fill(NaN, size(idxs))
    depths[.!isnothing.(idxs)] = alldepths[idxs[.!isnothing.(idxs)]]
    return depths
end
function getchanneldepths(session, probeid, X::LFPMatrix)
    channels = dims(X, Dim{:channel})|>collect
    getchanneldepths(session, probeid, channels)
end

function getunitdepths(session, probeid, units)
    metrics = getunitanalysismetrics(session)
    check = all(units .∈ (metrics.ecephys_unit_id,))
    if !check
        @error "Some units do not have corresponding metrics"
    end
    metrics = subset(metrics, :ecephys_unit_id, units)
    channels = metrics.ecephys_channel_id
    getchanneldepths(session, probeid, channels)
end

getdim(X::AbstractDimArray, dim) = dims(X, dim).val
gettimes(X::AbstractDimArray) = getdim(X, Ti)


function sortbydepth(session, probeid, LFP::AbstractDimArray)
    depths = getchanneldepths(session, probeid, getchannels(LFP))
    indices = Array{Any}([1:size(LFP, i) for i in 1:length(size(LFP))])
    indices[findfirst(isa.(dims(LFP), Dim{:channel}))] = sortperm(depths)
    return LFP[indices...]
end


function rectifytime(X::AbstractDimArray; tol=5) # tol gives significant figures for rounding
    times = gettimes(X)
    step = times |> diff |> mean
    err = times |> diff |> std
    if err > step/10.0^(-tol)
        @warn "Time step is not approximately constant, skipping rectification"
    else
        step = round(step; sigdigits=tol)
        t0, t1 = round.(extrema(times); sigdigits=tol+1)
        times = t0:step:t1+(10*step)
        times = times[1:size(X, Ti)] # Should be ok?
    end
    @assert length(times) == size(X, Ti)
    X = set(X, Ti => times)
end

function stimulusepochs(session, stim)
    stimtable = getepochs(session, stim)
    stimtable.interval = [a..b for (a, b) in zip(stimtable.start_time, stimtable.stop_time)]
    return stimtable
end
function stimulusintervals(session, stim)
    stimtable = getstimuli(session, stim)
    stimtable.interval = [a..b for (a, b) in zip(stimtable.start_time, stimtable.stop_time)]
    return stimtable
end

function gaborintervals(session)
    stimtable = getstimuli(session, "gabors")
    stimtable.combined_pos = sqrt.(Meta.parse.(stimtable.x_position).^2 .+ Meta.parse.(stimtable.y_position).^2) # Radial position of the gabor stimulus
    stimtable.interval = [a..b for (a, b) in zip(stimtable.start_time, stimtable.stop_time)]
    return stimtable
end

function radialgaborseries(session, times)
    stimtable = getstimuli(session, "gabors")
    stimtable.combined_pos = sqrt.(Meta.parse.(stimtable.x_position).^2 .+ Meta.parse.(stimtable.y_position).^2) # Radial position of the gabor stimulus
    gaborseries = DimArray(zeros(length(times)), (Ti(times),))
    for pos in unique(stimtable.combined_pos)
        gabortimes = [a..b for (a, b, x) in zip(stimtable.start_time, stimtable.stop_time, stimtable.combined_pos) if x == pos]
        for ts in gabortimes
            gaborseries[Ti(ts)] .= pos
        end
    end
    return gaborseries
end

function alignlfp(session, X, ::Val{:gabors}; x_position=nothing, y_position=nothing)
    gaborstim = gaborintervals(session)
    X = rectifytime(X)
    isnothing(x_position) || (gaborstim = gaborstim[Meta.parse.(gaborstim.x_position) .== x_position, :])
    isnothing(y_position) || (gaborstim = gaborstim[Meta.parse.(gaborstim.y_position) .== y_position, :])
    _X = [X[Ti(g)] for g in gaborstim.interval]
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X] # Catch any that are one sample too long
    # _X = DimArray(mean(collect.(_X)), (Ti(step(dims(X, Ti)):step(dims(X, Ti)):step(dims(X, Ti))*minimum(size.(_X, Ti))), dims(X, Dim{:channel})))
    return _X
end

function alignlfp(session, X, ::Val{:static_gratings})
    stim = stimulusintervals(session, "static_gratings")
    stim = stim[stim.start_time .> minimum(dims(X, Ti)), :]
    stim = stim[stim.stop_time .< maximum(dims(X, Ti)), :]
    X = rectifytime(X)
    _X = [X[Ti(g)] for g in stim.interval]
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X]
    return _X
end

"""
For flashes alignment, `trail=false` will return only the data from within the flash period. `trail=onset` will return the data from the onset of the flash to the onset of the flash through to the onset of the next flash. `trail=offset` will return the data from the offset of the flash to the onset of the next flash.
"""
function alignlfp(session, X, ::Val{:flashes}; trail=true)
    is = stimulusintervals(session, "flashes")
    if trail == :onset
        onsets = is.start_time
        is = [onsets[i]..onsets[i+1] for i in 1:length(onsets)-1]
    elseif trail == :offset
        offsets = is.stop_time
        onsets = is.start_time[2:end]
        is = [offsets[i]..onsets[i] for i in 1:length(offsets)-1]
    else
        is = is.interval
    end
    X = rectifytime(X)
    _X = [X[Ti(g)] for g in is]
    _X = [x[1:minimum(size.(_X, Ti)), :] for x in _X] # Catch any that are one sample too long
    return _X
end

# """
# For spontaneous alignment, we take each whole spontaneous interval as a trial
# """
# function alignlfp(session, X, ::Val{:spontaneous})
#     is =stimulusintervals(session, "spontaneous").interval
#     X = rectifytime(X)
#     _X = [X[Ti(g)] for g in is]
#     return _X
# end

alignlfp(session, X, stimulus::Union{String, Symbol}="gabors"; kwargs...) = alignlfp(session, X, stimulus|>Symbol|>Val; kwargs...)




bandpass(; kwargs...) = x -> bandpass(x; kwargs...)

function bandpass(x::LFPVector; pass, designmethod=Butterworth(4))
    t = dims(x, Ti)
    T = t isa AbstractRange ? step(t) : t |> collect |> diff |> mean
    fs = 1.0/T
    y = filtfilt(digitalfilter(Bandpass(pass...; fs), designmethod), x)
    return DimArray(y, dims(x))
end

function bandpass(X::LFPMatrix, dim=Dim{:channel}; kwargs...)
    Y = similar(X)
    Threads.@threads for x in dims(X, dim)
        Y[dim(At(x))] .= bandpass(X[dim(At(x))]; kwargs...)
    end
    return Y
end

thetafilter(args...; pass=[2, 8], kwargs...) = bandpass(args...; pass, kwargs...)
gammafilter(args...; pass=[30, 400], kwargs...) = bandpass(args...; pass, kwargs...)
broadbandfilter(args...; pass=[10, 400], kwargs...) = bandpass(args...; pass, kwargs...)
narrowbandfilter(args...; pass=[50, 60], kwargs...) = bandpass(args...; pass, kwargs...)


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
function thetafeature(x::AbstractVector, fs=nothing; τ=50, durprop=0.5)
    # ! fs not implemented. τ in time steps
    x = x[1:round(Int, length(x)*durprop)]
    return selfmutualinfo(x, [τ]) |> first
end


function thetafeature(x::LFPVector; kwargs...)
    fs = 1.0/step(dims(x, Ti))
    thetafeature(x, fs; kwargs...)
end



"""
Calculate a feature profile for each channel in each region
"""
function stimuluspartition(session, probeids, structures, stim; inbrain=200, times=nothing, kwargs...)
    Y = Vector{AbstractVector}([])
    stim == "spontaneous" && return spontaneouspartition(session, probeids, structures; inbrain=200, kwargs...)
    epochs = getepochs(session, stim)
    for p in eachindex(probeids)
        if isnothing(times)
            Z = []
            for e in axes(epochs, 1)
                times = epochs[e, :].start_time..epochs[e, :].stop_time
                X = getlfp(session, probeids[p], structures[p]; inbrain, times) |> rectifytime
                X = alignlfp(session, X, stim; kwargs...)
                push!(Z, X)
            end
            X = vcat(Z...)
        else
            X = getlfp(session, probeids[p], structures[p]; inbrain, times) |> rectifytime
            X = alignlfp(session, X, stim; kwargs...)
        end
        push!(Y, X)
    end
    return Y
end

function spontaneouspartition(session, probeids, structures; inbrain=200, mindur=30, kwargs...)
    times = stimulusintervals(session, "spontaneous").interval
    idxs = (diff.(extrema.(times).|>collect).|>first) .> mindur
    Y = Vector{AbstractVector}([])
    for p in eachindex(probeids)
        @info "Loading LFP for probe $p"
        X = []
        for (i, t) in enumerate(times)
            try
                _X = getlfp(session, probeids[p], structures[p]; inbrain, times=t) |> rectifytime
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


function spontaneouspartition(session, probeids, structures, duration; inbrain=200)
    epoch = getepochs(session, "spontaneous")[2, :]
    times = epoch.start_time..epoch.stop_time
    Y = Vector{AbstractVector}([])
    for p in eachindex(probeids)
        X = getlfp(session, probeids[p], structures[p]; inbrain, times) |> rectifytime
        _X = []
        t = minimum(dims(X, Ti))
        while t + duration < maximum(dims(X, Ti))
            push!(_X, X[Ti(t..(t+duration))])
            t += duration
        end
        push!(Y, _X)
    end
    return Y
end


"""
Calculate the thetafeature for each stimulus presentation
"""
function thetafeature(Y::Vector{AbstractVector}; kwargs...) # Input formatted as stimuluspartition()
    t = [[([mean(dims(c, Ti)) for c in eachcol(s)]) for s in y] for y in Y]
    F = [[([thetafeature(c; kwargs...) for c in eachcol(s)]) for s in y] for y in Y]
    return F, t
end

# function meanfeature(Y::AbstractVector{AbstractVector{AbstractVector}})
#     [[([mean(dims(c, Ti)) for c in eachcol(s)]) for s in y] for y in Y]
# end


function ica(X::LFPMatrix, k=ceil(Int, size(X, 2)/2))
    I = fit(ICA, collect(X)', k; maxiter=1000)
    _X = DimArray(predict(I, collect(X)')', (dims(X, Ti), Dim{:channel}(1:size(I, 2))))
end

function pca(X::LFPMatrix)
    P = fit(PCA, collect(X)')
    _X = DimArray(predict(P, collect(X)')', (dims(X, Ti), Dim{:channel}(1:size(P, 2))))
    v = principalvars(P)
    return _X, v
end

DSP.hilbert(X::LFPMatrix) = mapslices(hilbert, X, dims=Ti)
DSP.hilbert(X::LFPVector) = DimArray(hilbert(X|>Array), dims(X); refdims=refdims(X))

function waveletfreqs(t; moth=Morlet(2π), β=1, Q=32)
    n = length(t)
    fs = 1.0./step(t) # Assume rectified time dim
    W = ContinuousWavelets.computeWavelets(n, wavelet(moth; β, Q);)[1]
    freqs = getMeanFreq(W, fs)
    freqs[1] = 0
    return freqs
end

function _wavelettransform(x::AbstractVector; moth, β, Q) # β = 1 means linear in log space
    c = wavelet(moth; β, Q);
    res = ContinuousWavelets.cwt(x, c)
end

function _wavelettransform(x::LFPVector; rectify=true, moth=Morlet(2π), β=1, Q=32)
    rectify && (x = rectifytime(x))
    res = _wavelettransform(x|>Array; moth=Morlet(2π), β=1, Q=32)
    t = dims(x, Ti)
    freqs = waveletfreqs(t; moth, β, Q)
    return DimArray(res, (t, Dim{:frequency}(freqs)); metadata=DimensionalData.metadata(x), refdims=refdims(x))
end

function _wavelettransform(x::LFPVector, ::Val{:mmap}; window=50000, kwargs...)
    md = DimensionalData.metadata(x)
    rd = DimensionalData.refdims(x)
    x = rectifytime(x)
    𝓍 = _slidingwindow(x, window; tail=:overlap)
    t = dims(x, Ti)
    e = step(t)/2
    freqs = waveletfreqs(dims(𝓍[1], Ti); kwargs...)
    sz = (length(t), length(freqs))
    fname = tempname()
    s = open(fname, "w+"); write.((s,), sz)
    W = mmap(s, Matrix{ComplexF32}, sz)
    res = DimArray(W, (t, Dim{:frequency}(freqs)); metadata=(; md..., file=fname), refdims=rd)
    threadlog, threadmax = (0, length(𝓍))
    @withprogress name="Wavelet transform" begin
        for _x in 𝓍
            subres = _wavelettransform(_x; rectify=false, kwargs...)
            tx = extrema(dims(subres, 1))
            fx = extrema(dims(subres, 2))
            tilims = Interval{:closed, :closed}(tx[1]-e, tx[2]+e)
            flims = Interval{:closed, :closed}(fx[1]-e, fx[2]+e)
            res[Ti(tilims), Dim{:frequency}(flims)] .= subres
            if threadmax > 1
                Threads.threadid() == 1 && (threadlog += 1)%1 == 0 && @logprogress threadlog/threadmax
            end
        end
    end
    close(s)
    return res
end

_wavelettransform(x, s::Symbol; kwargs...) = _wavelettransform(x, Val(s); kwargs...)

wavelettransform(x::LFPVector, args...; kwargs...) = abs.(_wavelettransform(x, args...; kwargs...))


function fooof(p::LogWaveletMatrix, freqrange=[1.0, 300.0])
    ffreqs = 10.0.^collect(dims(p, Dim{:logfrequency}))
    freqrange = py"[$(freqrange[1]), $(freqrange[2])]"o
    spectrum = vec(collect(p))
    fm = PyFOOOF.FOOOF(peak_width_limits=py"[0.5, 50.0]"o, max_n_peaks=10, aperiodic_mode="knee", peak_threshold=0.5)
    # if fm.aperiodic_params_[2] < 0.0
    #     fm = PyFOOOF.FOOOF(peak_width_limits=py"[0.5, 20.0]"o, max_n_peaks=4, aperiodic_mode="fixed")
    # end
    # fm.report(freqs, spectrum, freqrange)
    fm.add_data(ffreqs, spectrum, freqrange)
    fm.fit()
    return fm
end

function logaperiodicfit(p::LogWaveletMatrix, freqrange=[1.0, 300.0], args...; doplot=false)
    fm = fooof(p, freqrange, args...)
    if doplot != false
        p = fm.plot(; plt_log=true, file_name=doplot, save_fig=true)
    end
    # * The aperiodic model, as described in doi.org/10.1038/s41593-020-00744-x
    b, k, χ = fm.aperiodic_params_
    k = max(k, 0.01)
    # b, k, χ = length(ps) == 3 ? ps : (ps[1], 0.0, ps[2])
    L = f -> 10.0.^(b - log10(k + (10.0^(f))^χ)) # Expects log frequency values
end

aperiodicfit(args...) = f -> (logaperiodicfit(args...)(log10(f)))



function _fooofedwavelet(res::LogWaveletMatrix, args...; kwargs...)
    psd = mean(res, dims=Ti)
    ffreqs = dims(psd, Dim{:logfrequency}) |> collect
    L = logaperiodicfit(psd, args...; kwargs...)
end
_fooofedwavelet(res::WaveletMatrix) = _fooofedwavelet(convert(LogWaveletMatrix, res))

function fooofedwavelet!(res::LogWaveletMatrix, freqrange=[1.0, 300.0]; kwargs...)
    psd = mean(res, dims=Ti)
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
    l[ffreqs .< fmin] = psd[ffreqs .< fmin] # Ensures no funky behaviour outside of the fit bounds
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

function aperiodicfit(psd::PSDVector, freqrange=[1.0, 300.0])
    ffreqs = dims(psd, Dim{:frequency}) |> collect
    freqrange = py"[$(freqrange[1]), $(freqrange[2])]"o
    spectrum = vec(collect(psd))
    fm = PyFOOOF.FOOOF(peak_width_limits=py"[0.5, 50.0]"o, max_n_peaks=10, aperiodic_mode="knee", peak_threshold=0.5)
    fm.add_data(ffreqs, spectrum, freqrange)
    fm.fit()
    b, k, χ = fm.aperiodic_params_
    k = max(k, 0.01)
    L = f -> 10.0.^(b - log10(k + (f)^χ))
end

function fooofedspectrum!(psd::PSDMatrix, freqrange=[1.0, 300.0]; kwargs...)
    L = aperiodicfit.(eachcol(psd), (freqrange,); kwargs...)
    ffreqs = dims(psd, Dim{:frequency}) |> collect
    L = hcat([l.(ffreqs) for l in L]...)
    psd .= psd .- L
end



function wavelettransform(X::LFPMatrix; kwargs...)
    res = []
    threadlog, threadmax = (0, size(X, 2)/Threads.nthreads())
    @withprogress name="Wavelet transform" begin
        Threads.@threads for c in 1:axes(X, 2)
            push!(res, wavelettransform(X[:, c]; kwargs...))
            Threads.threadid() == 1 && (threadlog += 1)%1 == 0 && @logprogress threadlog/threadmax
        end
    end
    ti = dims(res[1], 1)
    freq = dims(res[1], 2)
    channel = dims(X, 2)
    return DimArray(cat(collect.(res)..., dims=3), (ti, freq, channel))
end



mmapwavelettransform(x::LFPVector; kwargs...) = wavelettransform(x, :mmap; kwargs...)


function wavelettransform!(res::Dict, LFP::LFPMatrix; window=false, kwargs...)
    threadlog, threadmax = (0, size(LFP, 2)/Threads.nthreads())
    @withprogress name="Wavelet transform" begin
        for c in axes(LFP, 2)
            x = LFP[:, c]
            if window > 0
                # * Window the time series and calculate the wavelet transform in manageable chunks
                𝓍 = _slidingwindow(x, window)
                t = dims(x, Ti)
                freqs = waveletfreqs(dims(𝓍[1], Ti); kwargs...)
                sz = (length(t), length(freqs))
                s = open(tempname(), "w+"); write.((s,), sz)
                W = mmap(s, Matrix{Float32}, sz)
                _res = DimArray(W, (t, Dim{:frequency}(freqs)))
                for _x in 𝓍
                    subres = wavelettransform(_x; kwargs...)
                    tx = extrema(dims(subres, 1))
                    fx = extrema(dims(subres, 2))
                    tilims = Interval{:closed, :closed}(tx[1]-eps(), tx[2]+eps())
                    flims = Interval{:closed, :closed}(fx[1]-eps(), fx[2]+eps())
                    _res[Ti(tilims), Dim{:frequency}(flims)] .= subres
                end
            else
                push!(res, dims(LFP, 2)[c]=>wavelettransform(x; kwargs...)) # Doesnt actually write to file
            end
            Threads.threadid() == 1 && (threadlog += 1)%1 == 0 && @logprogress threadlog/threadmax
        end
    end
end


TimeseriesSurrogates.surrogate(x::LFPVector, S::Surrogate; kwargs...) = (y = deepcopy(x); y .= surrogate(x|>collect.|>Float64, S; kwargs...).|>eltype(y); y)

function TimeseriesSurrogates.surrogenerator(x::LFPVector, S::IAAFT)
    sg = surrogenerator(x|>collect.|>Float64, S)
    return ()->DimArray(sg().|>eltype(x), dims(x))
end
function TimeseriesSurrogates.surrogenerator(X::LFPMatrix, S::IAAFT)
    sg = [surrogenerator(x|>collect.|>Float64, S) for x in eachcol(X)]
    function out()
        Y = deepcopy(X)
        [y .= s().|>eltype(X) for (s, y) in zip(sg, eachcol(Y))]
        return Y
    end
    return out
end


function powerspectra(t, X; n=5000, window=DSP.Windows.hanning)
    Δt = t[2] - t[1]
    @assert all(Δt .≈ diff(t))
    fp = x -> welch_pgram(x, n; fs=1/Δt, window)
    P = [fp(Array(x)) for x ∈ eachcol(X)]
    𝑓 = P[1].freq
    psd = hcat([p.power for p ∈ P]...)
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
