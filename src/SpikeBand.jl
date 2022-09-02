using SparseArrays
using ProgressLogging
using Distances
using Clustering

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


function _getspikes(units, times, amplitudes, _times, bin, rectify)
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
    @withprogress name="catch22" begin
    for u in eachindex(units)
        _times = times[u]
        _amplitudes = amplitudes[u]
        for t in eachindex(_times) # Hella slow
            spikes[Ti(Near(_times[t])), Dim{:unit}(u)] += _amplitudes[t]
        end
        @logprogress u/length(units)
    end
    end
    return spikes
end

"""
We combine the spike times and spike amplitudes into one sparse array, using a given bin width.
"""
function getspikes(S, timebounds=nothing; bin=1e-4, rectify=4)
    times = S |> getspiketimes
    units = times |> keys |> collect
    times = times |> values |> collect
    amplitudes = S |> getspikeamplitudes
    amplitudes = getindex.((amplitudes,), units)
    @assert all(length.(times) .== length.(amplitudes))
    # Set up the sparse dim array, then allocate spikes to the nearest time index
    _times = isnothing(times) ? vcat(_times...) : timebounds
    _getspikes(units, times, amplitudes, _times, bin, rectify)
end

function getspikes(S, stimulus::String; kwargs...)
    timebounds = first(getstimulustimes(S, stimulus)) # The first flashes stimulus set
    getspikes(S, timebounds; kwargs...)
end

function getspikes(S, stimulus::String, structure::String; kwargs...)
    spikes = getspikes(S, stimulus; kwargs...)
    structures = getstructureacronyms(S, dims(spikes, Dim{:unit}))
    spikes = spikes[:, structures .== structure]
end

function getstructureacronyms(session::AbstractSession, units)
    unittable = getunitmetrics(session)
    acronyms = Vector{Any}(undef, size(units))
    [acronyms[i] = unittable[unittable.unit_id.==units[i], :ecephys_structure_acronym][1] for i ∈ 1:length(units)]
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
