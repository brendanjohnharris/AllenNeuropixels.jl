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



    include("./ReferenceAtlas.jl")
end