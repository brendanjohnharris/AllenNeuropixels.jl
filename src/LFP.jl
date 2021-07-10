using DimensionalData
using Parquet
using Tables
using PyCall


function timearray!(X::Parquet.Table) # There should be a time column, and everything else is a time series
    @assert :time ∈ keys(X)
    times = X.time
    arechannels = setdiff(1:length(keys(X)), [findfirst(keys(X).==:time)])
    channels = setdiff(keys(X), [:time])
    times = Ti(times)
    channels = Dim{:channel}(channels)
    X = DimArray(Tables.matrix(X)[:, arechannels], (times, channels))
end

# Always downsample 10x in time
function getdownsampledlfp(S::AbstractSession, probeid::Int)
    time=10
    space=1

    @assert(any(getprobeids(S) .== probeid), "Probe $probeid does not belong to session $(getid(S))")
    @assert(subset(getprobes(S), :id=>ByRow(==(probeid)))[!, :has_lfp_data][1], @error "Probe $probeid does not have LFP data")
    path = joinpath(datadir, "Ecephys", "session_"*string(getid(S)), "probe_"*string(probeid)*"_lfp.parquet");
    if !isfile(path)
        @info "Accessing LFP via the Allen SDK. This can take a few minutes."
        # pydata = S.pyObject.get_lfp(probeid)
        # We are going to save the bulky LFP data in parquet file, without all the metadata contained in the allen nwb file. We'll then load this into julia directly to save some time at the expense of duplicating data on the disk
        #py"""$(parquet).write_table($(pyarrow).Table.from_pandas($(S.pyObject).get_lfp($probeid).to_pandas()), $path)"""
        py"""$(parquet).write_table($(pyarrow).Table.from_pandas($(S.pyObject).get_lfp($probeid).to_pandas().iloc[::$time, ::$space]), $path)"""
    end
    @info "Loading from parquet file. Expect another minute"
    Float32.(timearray!(read_parquet(path)))
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


