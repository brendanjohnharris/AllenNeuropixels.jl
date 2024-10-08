module AllenNeuropixels
using PythonCall
using Preferences
using DataFrames
using Dates
using TimeZones
using CSV
using Requires
using IntervalSets
using Reexport
using PythonCall
using TimeseriesTools
@reexport using AllenNeuropixelsBase
import AllenNeuropixelsBase.rectifytime
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" @eval include("./Plots/Plots.jl")
    # @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" @eval using .Plots
end

include("./Windows.jl")
include("./LFP.jl")
include("./SpikeBand.jl")
include("./Bursts.jl")
include("./SpikeCluster.jl")

end
