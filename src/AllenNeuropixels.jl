module AllenNeuropixels
using PythonCall
using Preferences
using Conda
using DataFrames
using Dates
using TimeZones
using CSV
using PyCall
using Requires
using IntervalSets
using Reexport
ENV["PYTHON"] = "python3.9"
ENV["JULIA_PYTHONCALL_EXE"] = "python3.9"
using PythonCall
@reexport using AllenNeuropixelsBase
function __init__()
    run(`$(PyCall.python) -m pip install numpy`)
    run(`$(PyCall.python) -m pip install https://github.com/brendanjohnharris/AllenSDK/archive/refs/tags/v2.19.zip`)
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" @eval include("./Plots/Plots.jl")
    # @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" @eval using .Plots
end


include("./Windows.jl")
include("./LFP.jl")
include("./SpikeBand.jl")
include("./Bursts.jl")
include("./SpikeCluster.jl")

end
