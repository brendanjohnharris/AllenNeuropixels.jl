module AllenNeuropixels
using PyCall
using Preferences
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
    Conda.update() # You might need to delete the .julia/conda folder and rebuild python; allensdk has some tricky compatibility requirements.
    Conda.pip("install", "allensdk")
    copy!(allensdk, pyimport("allensdk"))
    copy!(brain_observatory, pyimport("allensdk.brain_observatory"))
    copy!(ecephys, pyimport("allensdk.brain_observatory.ecephys"))
    copy!(ecephys_project_cache, pyimport("allensdk.brain_observatory.ecephys.ecephys_project_cache"))
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=manifestpath)
end


function setdatadir(datadir::String)
    @set_preferences!("datadir" => datadir)
    @info("New default datadir set; restart your Julia session for this change to take effect")
end
const datadir = @load_preference("datadir", joinpath(pkgdir(AllenNeuropixels), "data/"))

function setmanifestpathr(manifestpath::String)
    @set_preferences!("manifestpath" => manifestpath)
    @info("New default manifestpath set; restart your Julia session for this change to take effect")
end
const manifestpath = @load_preference("manifestpath", joinpath(datadir, "manifest.json"))


function loaddataframe(file, dir=datadir)
    CSV.File(abspath(dir, file)) |> DataFrame
end
export convertdataframe


include("./EcephysCache.jl")


end