module AllenNeuropixels
using PyCall
using Conda
using DataFrames
using Dates
using TimeZones


const allensdk = PyNULL()
const brain_observatory = PyNULL()
const ecephys = PyNULL()
const ecephys_project_cache = PyNULL()
export allensdk, brain_observatory, ecephys, ecephys_project_cache

function __init__()
    Conda.pip_interop(true)
    #Conda.update()
    Conda.pip("install", "allensdk")
    copy!(allensdk, pyimport_conda("allensdk", "allendsk"))
    copy!(brain_observatory, pyimport_conda("allensdk.brain_observatory", "allendsk"))
    copy!(ecephys, pyimport_conda("allensdk.brain_observatory.ecephys", "allendsk"))
    copy!(ecephys_project_cache, pyimport_conda("allensdk.brain_observatory.ecephys.ecephys_project_cache", "allendsk"))
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=manifestdir)
end

const datadir = joinpath(pkgdir(AllenNeuropixels), "data/")
const manifestdir = joinpath(datadir, "manifest.json")


function convertdataframe(pdf, types)
    nc = 1:length(pdf.columns)
    nr = 1:length(pdf)
    columns = Symbol.(convert.((String,), getindex.((pdf.columns,), nc)))
    rows =  (convert.((Int,), getindex.((pdf.index,), nr)))
    df = DataFrame()
    df.id = rows
    M = pdf.values[:, :]
    for c ∈ nc
        v = Vector{types[c]}(undef, length(nr))
        println(c)
        for r ∈ nr
            if (try
                    isnan(convert(Float64, M[r, c]))
                catch
                    false
                end) && !(M[r, c] isa Float64)
                v = convert(Vector{Any}, v) # Prolly slow
                v[r] = convert(Float64, M[r, c])
            else
                types[c] <: AbstractVector ? v[r] = collect(M[r, c]) : v[r] = convert(types[c], M[r, c])
            end
        end
        df[!, columns[c]] = v
    end
    return df
end
export convertdataframe



include("./EcephysProject.jl")
include("./EcephysCache.jl")


end