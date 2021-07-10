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

    # WGLMakie works best from Julia in the terminal, not in VSCode
    #Page(exportable=true, offline=true)


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




    include("./ReferenceAtlas.jl")
    include("./LFP.jl")
end