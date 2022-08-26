using SparseArrays
using ProgressLogging

function downloadspikes(S::AbstractSession)
    _ = S.pyObject.spike_times
    _ = S.pyObject.spike_amplitudes
    return nothing
end





function getsessionpath(session::AbstractSession)
    path = joinpath(datadir, "Ecephys", "session_"*string(getid(session)), "session_"*string(getid(session))*".nwb");
end

getspiketimes(S::AbstractSession) = S.pyObject.spike_times
getspikeamplitudes(S::AbstractSession) = S.pyObject.spike_amplitudes


"""
We combine the spike times and spike amplitudes into one sparse array, using a given bin width.
"""
function getspikes(S::AbstractSession; bin=1e-4, rectify=4)
    times = S |> getspiketimes
    units = times |> keys |> collect
    times = times |> values |> collect
    amplitudes = S |> getspikeamplitudes
    amplitudes = getindex.((amplitudes,), units)
    @assert all(length.(times) .== length.(amplitudes))
    # Set up the sparse dim array, then allocate spikes to the nearest time index
    tmax = maximum(maximum.(times))
    tmin = minimum(minimum.(times))
    if rectify > 0
        tmax = round(tmax; digits=rectify)
        tmin = round(tmin; digits=rectify)
    end
    ts = tmin:bin:tmax
    spikes = spzeros(Float32, length(ts), length(units))
    spikes = goSparseDimArray(spikes, (Ti(ts), Dim{:unit}(units))) # SparseDimArray?
    @withprogress name="catch22" begin
    for u in eachindex(units)
        _times = times[u]
        _amplitudes = amplitudes[u]
        for t in eachindex(_times) # Hella slow
            spikes[Ti(Near(_times[t])), Dim{:unit}(u)] = _amplitudes[t]
        end
        @logprogress u/length(units)
    end
    end
    return spikes
end


function minspikediff(S::AbstractSession)
    diffs = S |> getspiketimes |> values .|> diff
    return minimum(minimum.(diffs)) # Just greater than 1e-4
end
