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
        @info "Slicetime Baybee"
        lfp = f["acquisition"][splitext(basename(path))[1]][splitext(basename(path))[1]*"_data"]["data"][channelidxs, timeidxs]
    elseif timeidxs isa Union{Int64, AbstractRange{Int64}}
        @info "Yeah we slice"
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

function getlfp(session, probeids::Vector{Int}, args...; kwargs...)
    LFP = [getlfp(session, probeid, args...; kwargs...) for probeid ∈ probeids]
end

function getchannels(data::DimArray)
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

getdim(X::AbstractDimArray, dim) = dims(X, dim).val
gettimes(X::AbstractDimArray) = getdim(X, Ti)
