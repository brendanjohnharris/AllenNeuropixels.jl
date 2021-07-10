using DimensionalData
using PyCall
using IntervalSets
using HDF5


function downloadlfp(S::AbstractSession, probeid::Int)
    @assert(any(getprobeids(S) .== probeid), "Probe $probeid does not belong to session $(getid(S))")
    @assert(subset(getprobes(S), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    _ = S.pyObject.get_lfp(probeid)
    return nothing
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
    f = h5open(getlfppath(session, probeid))
    channels = f["general"]["extracellular_ephys"]["electrodes"]["id"][:]
    close(f)
    return channels
end


"""
Get the lfp data for a probe,  providing *indices* for channels and times. See function below for indexing by channel ids and time values/intervals
"""
function _getlfp(session::AbstractSession, probeid::Int; channels=1:length(getlfpchannels(session, probeid)), times=1:length(getlfptimes(session, probeid)))
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
    timedata = timedata[times]

    channelids = getlfpchannels(session, probeid)
    channelids = channelids[channels]
    if (channels isa Union{Int64, AbstractRange{Int64}}) & (times isa Union{Int64, AbstractRange{Int64}}) # Can use HDF5 slicing
        @info "Slicetime Baybee"
        lfp = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][channels, times]
    else
        lfp = read(f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"])
        lfp = lfp[channels, times]
    end
    X = DimArray(lfp', (Ti(timedata),  Dim{:channel}(channelids)))
    close(f)
    return X
end


function isinvalidtime(session, probeids=getprobeids(session), times=NaN)
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
function getlfp(session::AbstractSession, probeid::Int; channels=getlfpchannels(session, probeid), times=OpenInterval(extrema(getlfptimes(session, probeid))...))
    channelids = channels
    timevals = times
    if isinvalidtime(session, probeid, timevals)
        @error "Requested LFP data contains an invalid time..."
    end

    channels = getlfpchannels(session, probeid)
    channels = indexin(channelids, channels)
    @assert length(channels) == length(channelids)

    times = getlfptimes(session, probeid)
    if !(timevals isa Interval) && length(timevals) == 2
        timevals = OpenInterval(timevals...)
    end
    if timevals isa Interval
        times = findall(times .∈ (timevals,))
    else
        times = indexin(timevals, times)
        @assert length(times) == length(timevals)
    end

    # See if we can convert to unitranges for faster HDF5 reading via slicing
    if collect(UnitRange(extrema(times)...)) == times
        times = UnitRange(extrema(times)...)
    end
    if collect(UnitRange(extrema(channels)...)) == channels
        channels = UnitRange(extrema(channels)...)
    end

    @info "Accessing LFP data"
    _getlfp(session, probeid; channels, times)
end

"""
If you want to downsample the LFP data, its quicker to use this function and then perform slicing afterwards (since getlfp() has to check all of the time coordinates you supply, which can be slow).
"""
function getdownsampledlfp(session, probeid; downsample=100, timerange=OpenInterval(extrema(getlfptimes(session, probeid))...))
    if !(timerange isa Interval) && length(timerange) == 2
        timerange = OpenInterval(timerange...)
    end
    timevals = getlfptimes(session, probeid)
    tidxs = timevals .∈ (timerange,)
    times = findfirst(tidxs):downsample:findlast(tidxs)
    _getlfp(session, probeid; times)
end