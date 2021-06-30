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
const mouse_connectivity_cache = PyNULL()
const ontologies_api = PyNULL()
const reference_space_cache = PyNULL()
export allensdk, brain_observatory, ecephys, ecephys_project_cache, mouse_connectivity_cache, ontolofies_api, reference_space_cache

function __init__()
    Conda.pip_interop(true)
    Conda.update() # You might need to delete the .julia/conda folder and rebuild python; allensdk has some tricky compatibility requirements.
    Conda.pip("install", "allensdk")
    copy!(allensdk, pyimport("allensdk"))
    copy!(brain_observatory, pyimport("allensdk.brain_observatory"))
    copy!(ecephys, pyimport("allensdk.brain_observatory.ecephys"))
    copy!(ecephys_project_cache, pyimport("allensdk.brain_observatory.ecephys.ecephys_project_cache"))
    copy!(mouse_connectivity_cache, pyimport("allensdk.core.mouse_connectivity_cache"))
    copy!(ontologies_api, pyimport("allensdk.api.queries.ontologies_api"))
    copy!(reference_space_cache, pyimport("allensdk.core.reference_space_cache"))
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=ecephysmanifest)
end


function setdatadir(datadir::String)
    @set_preferences!("datadir" => datadir)
    @info("New default datadir set; restart your Julia session for this change to take effect")
end
const datadir = replace(@load_preference("datadir", joinpath(pkgdir(AllenNeuropixels), "data/")), "\\"=>"/")
const ecephysmanifest = replace(joinpath(datadir, "Ecephys", "manifest.json"), "\\"=>"/")
const mouseconnectivitymanifest = replace(joinpath(datadir, "MouseConnectivity", "manifest.json"), "\\"=>"/")
const referencespacemanifest = replace(joinpath(datadir, "ReferenceSpace", "manifest.json"), "\\"=>"/")

function loaddataframe(file, dir=datadir)
    CSV.File(abspath(dir, file)) |> DataFrame
end
export convertdataframe


include("./EcephysCache.jl")
include("./LFP.jl")
include("./MouseConnectivityCache.jl")
include("./Plotting.jl")
include("./Ontologies.jl")
include("./ReferenceSpace.jl")


end