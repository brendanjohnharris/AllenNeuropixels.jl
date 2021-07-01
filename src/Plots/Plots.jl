module Plots
    using WGLMakie
    using JSServe
    using Colors
    using Meshing
    using MeshIO
    using GeometryBasics
    using DataFrames
    import ..AllenNeuropixels as AN

    # This works best from Julia in the terminal, not in VSCode
    Page(exportable=true, offline=true)



    include("./ReferenceAtlas.jl")
end