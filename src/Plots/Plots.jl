module Plots
    using ..Makie
    using JSServe
    using Colors
    using Meshing
    using MeshIO
    using Markdown
    using GeometryBasics
    using DataFrames
    using FileIO
    import ..AllenNeuropixels as AN

    function htmlexport(fig, file::String="plot.html")
        open(file, "w") do io
                println(io, """
                <html>
                    <head>
                    </head>
                    <body>
                """)
                show(io, MIME"text/html"(), Page(exportable=true, offline=true))
                show(io, MIME"text/html"(), fig)
                println(io, """
                    </body>
                </html>
                """)
        end
    end
    exporthtml = htmlexport

    Makie.convert_arguments(x::AN.LFPVector) = (dims(x, Ti)|>collect, x|>Array)
    Makie.convert_single_argument(x::AN.LFPVector) = (dims(x, Ti)|>collect, x|>Array)

    GeometryBasics.decompose(x::AN.DimensionalData.AbstractDimArray) = ((dims(x).|>collect)..., x.data)

    dimname(x::AN.DimensionalData.AbstractDimArray, dim) = dims(x, dim)|>name|>string

    formataxes(x::AN.DimensionalData.AbstractDimArray{T, 2} where T) = (xlabel=dimname(x, 1), ylabel=dimname(x, 2))
    formataxes(x::AN.DimensionalData.AbstractDimArray{T, 1} where T) = (xlabel=dimname(x, 1),)

    include("./ReferenceAtlas.jl")
    include("./LFP.jl")
    include("SpikeBand.jl")
end
