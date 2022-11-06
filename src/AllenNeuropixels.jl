module AllenNeuropixels
using PyCall
using Preferences
using Conda
using DataFrames
using Dates
using TimeZones
using CSV
using Requires
using IntervalSets

# ! For the best chance at a working install:
# 1. Manually install allensdk in python3.7
# ENV["PYTHON"]="/usr/bin/python3.7"

const pynwb = PyNULL()
const allensdk = PyNULL()
const brain_observatory = PyNULL()
const ecephys = PyNULL()
const stimulusmapping = PyNULL()
const ecephys_project_cache = PyNULL()
const ecephys_project_api = PyNULL()
const ephys_features = PyNULL()
const brain_observatory_cache = PyNULL()
const stimulus_info = PyNULL()
const mouse_connectivity_cache = PyNULL()
const ontologies_api = PyNULL()
const reference_space_cache = PyNULL()
const reference_space = PyNULL()
export allensdk, brain_observatory, ecephys, ecephys_project_cache, mouse_connectivity_cache, ontologies_api, reference_space_cache, reference_space

function __init__()
    # Conda.pip_interop(true)
    # Conda.update() # You might need to delete the .julia/conda folder and rebuild PyCall; allensdk has some tricky compatibility requirements.
    # Conda.pip("install", "allensdk")

    copy!(allensdk, pyimport("pynwb"))
    copy!(allensdk, pyimport("allensdk"))
    copy!(brain_observatory, pyimport("allensdk.brain_observatory"))
    copy!(stimulus_info, pyimport("allensdk.brain_observatory.stimulus_info"))
    copy!(ecephys, pyimport("allensdk.brain_observatory.ecephys"))
    copy!(stimulusmapping, pyimport("allensdk.brain_observatory.ecephys.stimulus_analysis.receptive_field_mapping"))
    copy!(ecephys_project_cache, pyimport("allensdk.brain_observatory.ecephys.ecephys_project_cache"))
    copy!(ecephys_project_api, pyimport("allensdk.brain_observatory.ecephys.ecephys_project_api"))
    copy!(ephys_features, pyimport("allensdk.ephys.ephys_features"))
    copy!(brain_observatory_cache, pyimport("allensdk.core.brain_observatory_cache"))
    copy!(mouse_connectivity_cache, pyimport("allensdk.core.mouse_connectivity_cache"))
    copy!(ontologies_api, pyimport("allensdk.api.queries.ontologies_api"))
    copy!(reference_space_cache, pyimport("allensdk.core.reference_space_cache"))
    copy!(reference_space, pyimport("allensdk.core.reference_space"))

    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=ecephysmanifest)

    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" @eval include("./Plots/Plots.jl")
    # @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" @eval using .Plots
end


# function setpython37(python39::String)
#     @set_preferences!("python39" => python39)
#     @info("New default python executable set; restart your Julia session for this change to take effect")
# end
# const python39 = @load_preference("python39", "")

function setdatadir(datadir::String)
    @set_preferences!("datadir" => datadir)
    @info("New default datadir set; restart your Julia session for this change to take effect")
end
const datadir = replace(@load_preference("datadir", joinpath(pkgdir(AllenNeuropixels), "data/")), "\\"=>"/")
const ecephysmanifest = replace(joinpath(datadir, "Ecephys", "manifest.json"), "\\"=>"/")
const brainobservatorymanifest = replace(joinpath(datadir, "BrainObservatory", "manifest.json"), "\\"=>"/")
const mouseconnectivitymanifest = replace(joinpath(datadir, "MouseConnectivity", "manifest.json"), "\\"=>"/")
const referencespacemanifest = replace(joinpath(datadir, "ReferenceSpace", "manifest.json"), "\\"=>"/")

function loaddataframe(file, dir=datadir)
    CSV.File(abspath(dir, file)) |> DataFrame
end
export convertdataframe


include("./EcephysCache.jl")
include("./BrainObservatory.jl")
include("./Streaming.jl")
include("./LFP.jl")
include("./SparseDimArray.jl")
include("./SpikeBand.jl")
include("./GammaBursts.jl")
include("./SpikeCluster.jl")
include("./MouseConnectivityCache.jl")
include("./Ontologies.jl")
include("./ReferenceSpace.jl")


end
