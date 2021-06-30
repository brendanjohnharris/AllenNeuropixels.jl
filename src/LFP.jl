using DimensionalData

function getlfp(S::AbstractSession, probeid::Int)
    @assert(any(getprobeids(S) .== probeid), "Probe $probeid does not belong to session $(getid(S))")
    @assert(subset(getprobes(S), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    pydata = S.pyObject.get_lfp(probeid)
    data = pydata.data
    times = Ti(pydata.indexes.get("time").data)
    channels = Dim{:channel}(pydata.indexes.get("channel").data)
    A = DimArray(data, (times, channels))
    #! Do the download, but read straight from nwb file
end


function getlfp(S::AbstractSession)
    lfp_data = 1
    for i ∈ getprobeids(S)
        # this_probe_id = session_7325.probes.index[i]
        # this_lfp = session.get_lfp(this_probe_id) # start downloading lfp data in␣
        # ,!the same data directory
        # lfp_data.append(this_lfp)
    end
end