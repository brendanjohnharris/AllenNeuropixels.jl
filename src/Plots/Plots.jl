module Plots
    using Makie
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

    include("./ReferenceAtlas.jl")
    include("./LFP.jl")
    include("SpikeBand.jl")
end
