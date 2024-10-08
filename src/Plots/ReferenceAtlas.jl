
rotatereferenceatlas = x -> reverse(permutedims(x, (1, 3, 2)), dims = (3,))
# Plot the reference volume and probe locations
function plotreferencevolume(S; dotext = true, dostructures = true, ids = :targets,
                             size = (1920, 1080), shading = Makie.automatic, kwargs...)
    channels = AN.getchannels(S)
    vol, info = AN.gettemplatevolume()
    vol = Array{Float16}(vol)
    grad = cgrad(cgrad([RGBA(1.0, 1.0, 1.0, 0.0), RGBA(0.0, 0.0, 0.0, 1.0)]), 100)
    #coords = [(1:s).*25 for s ∈ size(vol)] # 25 μm size
    #vol = reverse(permutedims(vol, [1 3 2]), dims=(1,3))
    #coords = [1:x for x ∈ size(vol)].*50
    s = Scene(backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.0), size = size,
              transparency = true)
    coords = [1:s for s in size(rotatereferenceatlas(vol)) ./ 2]
    coords[3] = .-coords[3]

    markercolor = :gray
    if ids isa Symbol
        if ids == :all
            # Get and plot a bunch of structures
            #ids = AN.getallstructureids()
            anns = unique(["FRP", "MO", "SS", "GU", "VISC", "AUD", "ACA", "PL", "ILA",
                              "ORB", "AI", "RSP", "PTLp", "TEa", "PERI", "ECT", "OLF",
                              "VISp", "VISl", "VISrl", "VISam", "VISpm", "VIS", "VISal",
                              "VISmma", "VISmmp", "VISli", "LGd", "LD", "LP", "VPM", "TH",
                              "MGm", "MGv", "MGd", "PO", "LGv", "VL",
                              "VPL", "POL", "Eth", "PoT", "PP", "PIL", "IntG", "IGL", "SGN",
                              "VPL", "PF", "RT", "CA1", "CA2", "CA3", "DG", "SUB", "POST",
                              "PRE", "ProS", "HPF", "MB", "SCig", "SCiw", "SCsg", "SCzo",
                              "PPT", "APN", "NOT", "MRN", "OP", "LT", "RPF", "CP"])
        elseif ids == :corticaltargets
            anns = unique(["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"])
        elseif ids == :targets # As listed in the visual coding white paper
            anns = unique(["VISp",
                              "VISl",
                              "VISrl",
                              "VISal",
                              "VISpm",
                              "VISam",
                              "CA1",
                              "CA3",
                              "DG",
                              "SUB",
                              "ProS",
                              "LGd",
                              "LP",
                              "APN"])
        else
            @error "`ids` keyword argument not valid "
        end
        t = AN.getstructuretree()
        D = t.get_id_acronym_map()
        ids = getindex.((D,), anns)
    end

    f = volume!(s, coords..., rotatereferenceatlas(vol)[1:2:end, 1:2:end, 1:2:end];
                algorithm = :mip, colorrange = extrema(vol), colormap = collect(grad),
                ticks = nothing, fxaa = true, transparency = true, kwargs...)
    #s.plots[1].attributes[:fxaa] = true

    chrome = FileIO.load(download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/E6BF3C_5A4719_977726_FCFC82.png"))
    for probeid in unique(channels.probe_id)
        (x, y, z) = AN.getprobecoordinates(S, probeid) ./ 50
        if !any(length.((x, y, z)) .== 0)
            f = meshscatter!(s, x, z, -y; markersize = 1.0, fxaa = true,
                             color = markercolor, matcap = chrome,
                             shading,
                             transparency = true, kwargs...)

            if dotext
                _, idx = findmax(z) # The highest unit
                text!(s, string(probeid), position = Vec3f0(x[idx], z[idx], -y[idx]),
                      space = :data, fontsize = 5, align = (:right, :bottom),
                      rotation = Quaternion((-0.3390201, 0.33899, 0.6205722, -0.620517)),
                      transparency = true)
            end
        end
    end
    #scale!(s, 1, 1, -1)
    if dostructures
        channels = subset(AN.getchannels(), :ecephys_probe_id => ByRow(!ismissing),
                          :ecephys_structure_id => ByRow(!ismissing),
                          :ecephys_structure_acronym => ByRow(x -> x !== "grey"))
        if isnothing(ids)
            ids = unique(channels[!, :ecephys_structure_id])
            #ids = ids[AN.getstructuretreedepth.(ids) .< 9]
        end
        for id in ids
            try
                try
                    mask, _ = AN.getstructuremask(id)
                catch
                    mask = AN.buildstructuremask(id)
                end
                mask = Array{Float64}(rotatereferenceatlas(mask)[1:2:end, 1:2:end, 1:2:end])
                # mc_algo = NaiveSurfaceNets(iso = 0, insidepositive = true)
                mc_algo = MarchingCubes(iso = 0, insidepositive = true)
                m = GeometryBasics.Mesh(mask, mc_algo;
                                        origin = [min(coords[i]...)
                                                  for i in 1:length(coords)],
                                        widths = [abs(-(extrema(coords[i])...))
                                                  for i in 1:length(coords)])
                c = AN.getstructurecolor(id)
                f = mesh!(s, m; color = RGBA(c.r, c.g, c.b, 0.41), fxaa = true,
                          shading, transparency = true, kwargs...)
            catch y
                @warn y
            end
        end
    end
    return s
end

function exportreferencevolume(S, file::String = "plot.html"; ids = :targets, kwargs...)
    open(file, "w") do io
        println(io,
                """
    <html>
        <head>
        </head>
        <body>
    """)
        show(io, MIME"text/html"(), Page(exportable = true, offline = true))
        fig = plotreferencevolume(S; dostructures = true, ids, kwargs...)
        rotate_cam!(fig, Vec3f0(0, 2.35, 0))
        show(io, MIME"text/html"(), fig)
        println(io,
                """
        </body>
    </html>
    """)
    end
end

function formattedreferencevolume(S, file::String = "plot.html";
                                  shading = Makie.automatic)
    exportreferencevolume(S, file;
                          dostructures = true,
                          ids = :targets,
                          show_axis = false,
                          shading, size = (1080, 1080))
end

# ? ---------------------------- # New functions --------------------------- ? #
function mesh2vertfaces(m)
    verts = hcat(collect.(m.position)...)' |> collect
    faces = [getproperty.(f, :i) .|> Int for f in GeometryBasics.faces(m)]
    faces = hcat(collect.(faces)...)' |> collect
    return (verts, faces)
end

ccftransform(x; scale = 1) = scale .* [1 0 0;
                                       0 0 1
                                       0 -1 0] * collect(x)
function ccftransform(m::GeometryBasics.Mesh; kwargs...)
    m.position .= ccftransform.(m.position; kwargs...)
    return m
end

"""
    plotbrainstructure!(ax, id; hemisphere)

Plot a brain structure with the given ID on the given axis.

# Arguments
- `ax::AbstractPlotting.Axis`: The axis to plot on.
- `id::Int`: The ID of the brain structure to plot.
- `hemisphere::Symbol`: The hemisphere to plot the structure on. Default is `:both`.

# Example
```
D = AllenNeuropixels.getstructureidmap()
id = D["root"] # The whole brain

f = Figure(; size = (1920, 1080))
ax = Axis3(f[1, 1]; aspect = :data)
p = AN.Plots.plotbrainstructure!(ax, id; hemisphere=:both)
```
"""
function plotbrainstructure!(ax, id; hemisphere = :both, alpha = 0.41, color = :atlas,
                             scale = 1,
                             kwargs...)
    m = ccftransform(ANB.getstructuremesh(id; hemisphere); scale)
    if color === :atlas
        c = AN.getstructurecolor(id)
    elseif color isa Symbol
        c = Makie.RGB(Makie.Colors.color_names[string(color)] ./ 255...)
    else
        c = color
    end
    mesh!(ax, m; color = RGBA(c.r, c.g, c.b, alpha), kwargs...)
end

"""
    plotbrain!(ax, S::AN.AbstractSession; dotext = :cortex, dostructures = true,
                ids = :targets, probestyle = :lines, dark = false, meshparams = ())

Plot a 3D representation of the brain with probes and structures.

# Arguments
- `ax::AbstractPlotting.Axis`: The axis to plot on.
- `S::AN.AbstractSession`: The session object containing the data.
- `dotext::Symbol=:cortex`: Whether to display text labels for the probes or the cortical structures.
- `dostructures::Bool=true`: Whether to plot the brain structures.
- `ids::Symbol=:targets`: Which structures to plot. Can be `:all`, `:corticaltargets`, `:targets`, or a vector of structure IDs.
- `probestyle::Symbol=:lines`: Whether to plot the probes as lines or meshscatter.
- `dark::Bool=false`: Whether to use a dark theme.
- `meshparams::NamedTuple=()`: Additional parameters to pass to the `mesh!` function.

# Returns
- A tuple containing the color observables and the probe plots.

# Example
```
S = ANB.VisualBehavior.Session(1067588044)
f = Figure(; size = (1920, 1080))
ax = Axis3(f[1, 1]; aspect = :data)
c, p = AN.Plots.plotbrain!(ax, S; dark = false)
```
"""
function plotbrain!(ax, S::AN.AbstractSession; dotext = :cortex, dostructures = true,
                    ids = :targets, probestyle = :meshscatter, dark = false,
                    meshparams = (), matcap = false, channels = nothing, markersize = 100.0,
                    fontsize = 25.0, shading = Makie.automatic, structurecolors = :atlas,
                    markercolor = :gold, hemisphere = :right)
    if !isnothing(channels)
        findchannels = channels
    else
        findchannels = nothing
    end
    channels = AN.getchannels(S)
    if !isnothing(findchannels)
        function findstructure(c)
            idxs = indexin(c, channels.id)
            channels[idxs, :].probe_id |> unique |>
            first
        end
        findchannels = Dict(findstructure.(findchannels) .=> findchannels)
    else
        findchannels = unique(channels.probe_id)
        findchannels = Dict(findchannels .=> [channels[channels.probe_id .== p, :].id
                             for p in findchannels])
    end
    vol, info = ANB.gettemplatevolume()

    if dark
        meshparams = (; fxaa = true,
                      shading,
                      transparency = true,
                      ssao = true,
                      #   interpolate = true,
                      #   ambient = Vec3f(0.1),
                      #   specular = Vec3f(2.0),
                      #   diffuse = Vec3f(0.1),
                      #   shininess = Float32(100.0),
                      meshparams...)
        rootalpha = 0.05
    else
        meshparams = (; fxaa = true,
                      shading,
                      transparency = true,
                      ssao = true,
                      #   interpolate = true,
                      #   ambient = Vec3f(0.5),
                      #   specular = Vec3f(0.2),
                      #   diffuse = Vec3f(0.2),
                      #   shininess = Float32(50.0),
                      meshparams...)
        rootalpha = 0.1
    end
    vol = Array{Float16}(vol)
    grad = cgrad(cgrad([RGBA(1.0, 1.0, 1.0, 0.0), RGBA(0.0, 0.0, 0.0, 1.0)]), 100)

    if ids isa Symbol
        if ids == :all
            # Get and plot a bunch of structures
            #ids = AN.getallstructureids()
            anns = unique(["FRP", "MO", "SS", "GU", "VISC", "AUD", "ACA", "PL", "ILA",
                              "ORB",
                              "AI", "RSP", "PTLp", "TEa", "PERI", "ECT", "OLF", "VISp",
                              "VISl",
                              "VISrl", "VISam", "VISpm", "VIS", "VISal", "VISmma", "VISmmp",
                              "VISli", "LGd", "LD", "LP", "VPM", "TH", "MGm", "MGv", "MGd",
                              "PO", "LGv", "VL",
                              "VPL", "POL", "Eth", "PoT", "PP", "PIL", "IntG", "IGL", "SGN",
                              "VPL", "PF", "RT", "CA1", "CA2", "CA3", "DG", "SUB", "POST",
                              "PRE", "ProS", "HPF", "MB", "SCig", "SCiw", "SCsg", "SCzo",
                              "PPT",
                              "APN", "NOT", "MRN", "OP", "LT", "RPF", "CP"])
        elseif ids == :corticaltargets
            anns = unique(["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"])
        elseif ids == :targets # As listed in the visual coding white paper
            anns = unique(["VISp",
                              "VISl",
                              "VISrl",
                              "VISal",
                              "VISpm",
                              "VISam",
                              "CA1",
                              "CA3",
                              "DG",
                              "SUB",
                              "ProS",
                              "LGd",
                              "LP",
                              "APN"])
        else
            @error "`ids` keyword argument not valid "
        end
    end
    if ids isa AbstractVector{<:String}
        anns = ids
    end
    D = ANB.getstructureidmap()
    ids = getindex.((D,), anns)
    # Plot the whole brain
    id = plotbrainstructure!(ax, 997; hemisphere = :both, meshparams..., alpha = rootalpha)

    if dostructures
        _channels = subset(AN.getchannels(), :ecephys_probe_id => ByRow(!ismissing),
                           :ecephys_structure_id => ByRow(!ismissing),
                           :ecephys_structure_acronym => ByRow(x -> !ismissing(x) &&
                                                                        x != "grey" &&
                                                                        x != "root"))
        if isnothing(ids)
            ids = unique(_channels[!, :ecephys_structure_id])
            #ids = ids[AN.getstructuretreedepth.(ids) .< 9]
        end
        if structurecolors isa Symbol || length(structurecolors) == 1
            structurecolors = fill(structurecolors, length(anns))
        end
        for (i, id) in (collect(enumerate(ids)))
            try
                plotbrainstructure!(ax, id; hemisphere, color = structurecolors[i],
                                    meshparams...)
            catch y
                @warn y
            end
        end
    end

    cp = map(enumerate(collect(findchannels))) do (i, fc)
        probeid = first(fc)
        xyz = AN.getchannelcoordinates(S, probeid)
        xyz = hcat(collect.(getindex.([xyz], last(fc)))...)
        (x, y, z) = eachrow(ccftransform(xyz))
        _c = nothing
        if !any(length.((x, y, z)) .== 0)
            if probestyle === :meshscatter && matcap == true
                chrome = FileIO.load(download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/E6BF3C_5A4719_977726_FCFC82.png"))
                # chrome = FileIO.load(download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/3B3C3F_DAD9D5_929290_ABACA8.png"))
                _c = Observable{Any}(chrome)
                _p = meshscatter!(ax, x, y, z; meshparams..., markersize,
                                  color = markercolor,
                                  matcap = _c, colormap = :bone)
            elseif probestyle === :meshscatter && matcap isa Matrix
                _c = Observable{Any}(matcap)
                _p = meshscatter!(ax, x, y, z; meshparams..., markersize,
                                  color = markercolor, matcap = _c, colormap = :bone)
            elseif probestyle === :meshscatter
                if markercolor isa Symbol || length(markercolor) == 1
                    _c = Observable{Any}(fill(markercolor, length(x)))
                elseif length(markercolor) == length(findchannels)
                    _c = Observable{Any}(fill(markercolor[i], length(x)))
                end
                _p = meshscatter!(ax, x, y, z; meshparams..., markersize,
                                  color = _c, colormap = :bone)
            elseif probestyle === :lines
                _c = Observable{Any}(fill(markercolor, length(x)))
                _p = lines!(ax, x, y, z; color = _c, linewidth = 25.0, colormap = :bone)
            end

            if dotext === :id
                _, idx = findmax(y) # The highest unit
                text!(ax, string(probeid) * "   ";
                      position = Vec3f(x[idx], y[idx], z[idx]),
                      space = :data, fontsize, align = (:right, :center),
                      rotation = Quaternion((-0.3390201, 0.33899, 0.6205722, -0.620517)),
                      transparency = true, overdraw = true)
            elseif dotext === :cortex
                structures = AN.getprobestructures(S)
                structures = Dict(p => s[contains.(s, ["VIS"])][1] for (p, s) in structures)
                _, idx = findmax(z) # The highest unit
                text!(ax, structures[probeid];
                      position = Vec3f(x[idx], y[idx], z[idx] + 250),
                      space = :data, fontsize, align = (:center, :center),
                      transparency = true, overdraw = true)
            end
        else
            _p = nothing
        end
        return (_c, _p)
    end

    if ax isa Axis
        ax.azimuth[] = 2.25
        ax.elevation[] = 0.2
    end
    return (first.(collect(findchannels)), first.(cp), last.(cp)) # The first are the color observables, the rest are the probes plots
end

function plotbrain!(ax, s::Integer, args...; kawrgs...)
    plotbrain!(ANB.VisualBehavior.Session(s), args...; kwargs...)
end

function label_connected_components(x::Vector{Int})
    y = deepcopy(x)
    v = first(x)
    y[1] = 1
    for i in eachindex(y)[2:end]
        if x[i] != v
            y[i] = y[i - 1] + 1
            v = x[i]
        else
            y[i] = y[i - 1]
        end
    end

    return y
end

function _layerplot(session::AN.AbstractSession, channels::AbstractVector{<:Int})
    channels = AN.sortbydepth(session, channels; method = :probe)
    depths = AN.getchanneldepths(session, channels; method = :probe)
    depths = Depth(depths)
    depths = first(rectify(depths))

    layerids, layers = AN.getchannellayers(session, channels) # By default, all channels along a probe. Subcortical structures don't get layers.

    unids = Int.(indexin(layers, unique(layers)))
    unids = label_connected_components(unids)

    dr = [0; diff(unids)] .> 0
    cs = [mean(depths[unids .== x]) for x in unique(unids)]
    return layers, depths, unids, dr, cs
end

function plotlayers!(ax, layers, depths)
    unids = Int.(indexin(layers, unique(layers)))
    unids = label_connected_components(unids)

    dr = [0; diff(unids)] .> 0
    cs = [mean(depths[unids .== x]) for x in unique(unids)]
    plotlayers!(ax, layers, depths, unids, dr, cs)
end

function plotlayers!(ax, layers, depths, unids, dr, cs)
    hlines!(ax, collect(depths[dr]) .- 0.5 * step(depths),
            color = (:white, 0.3))

    labels = [first(layers[unids .== x]) for x in unique(unids)]
    ax.yticks = (cs, labels)
end

function plotlayers!(ax, session::AN.AbstractSession, channels::AbstractVector{<:Int})
    layers, depths, unids, dr, cs = _layerplot(session, channels)
    plotlayers!(ax, layers, depths, unids, dr, cs)
end

function plotlayers!(ax, session::AN.AbstractSession, probeid::Integer)
    channels = AN.getlfpchannels(session, probeid)
    plotlayers!(ax, session, channels)
end
