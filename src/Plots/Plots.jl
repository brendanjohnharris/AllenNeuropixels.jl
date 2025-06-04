module Plots
using ..Makie
using TimeseriesTools
# using JSServe
using Colors
using Meshing
using MeshIO
using Markdown
using GeometryBasics
using DataFrames
using FileIO
import ..AllenNeuropixels as AN
import ..AllenNeuropixels: Chan, Unit, Depth, Log𝑓

# function htmlexport(fig, file::String="plot.html")
#     open(file, "w") do io
#             println(io, """
#             <html>
#                 <head>
#                 </head>
#                 <body>
#             """)
#             show(io, MIME"text/html"(), Page(exportable=true, offline=true))
#             show(io, MIME"text/html"(), fig)
#             println(io, """
#                 </body>
#             </html>
#             """)
#     end
# end
# exporthtml = htmlexport

Makie.convert_arguments(x::AN.LFPVector) = (dims(x, 𝑡) |> collect, x |> Array)
Makie.convert_single_argument(x::AN.LFPVector) = (dims(x, 𝑡) |> collect, x |> Array)

# GeometryBasics.decompose(x::AN.TimeseriesTools.AbstractToolsArray) = ((dims(x).|>collect)..., x.data)

# GeometryBasics.decompose(x::AN.TimeseriesTools.AbstractToolsArray, dims...) = (getindex.((dims(x).|>collect), dims)..., x.data[dims...])

dimname(x::AN.TimeseriesTools.AbstractToolsArray, dim) = dims(x, dim) |> name |> string

function formataxes(x::AN.TimeseriesTools.AbstractToolsArray{T, 2} where {T})
    (xlabel = dimname(x, 1), ylabel = dimname(x, 2))
end
function formataxes(x::AN.TimeseriesTools.AbstractToolsArray{T, 1} where {T})
    (xlabel = dimname(x, 1),)
end

include("./ReferenceAtlas.jl")
include("./LFP.jl")
include("SpikeBand.jl")
end
