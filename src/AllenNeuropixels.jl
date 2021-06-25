module AllenNeuropixels
using PyCall
using Conda
using DataFrames
using Dates
using TimeZones
using CSV


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
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=manifestpath)
end

const datadir = joinpath(pkgdir(AllenNeuropixels), "data/")
const manifestpath = joinpath(datadir, "manifest.json")


function loaddataframe(file, dir=datadir)
    df = CSV.File(abspath(dir, file)) |> DataFrame
end
export convertdataframe



include("./EcephysProject.jl")
include("./EcephysCache.jl")


end