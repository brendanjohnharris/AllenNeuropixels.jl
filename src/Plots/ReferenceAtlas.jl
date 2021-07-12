

rotatereferenceatlas = x -> reverse(permutedims(x, (1, 3, 2)), dims=(3,))
# Plot the reference volume and probe locations
function plotreferencevolume(S; dotext=true, dostructures = true, ids=:targets, resolution=(1920, 1080), kwargs...)
    #! Check probe locations are correct
    channels = AN.getchannels(S)
    vol, info = AN.gettemplatevolume()
    vol = Array{Float16}(vol)
    grad = cgrad(cgrad([RGBA(1.0, 1.0, 1.0, 0.0), RGBA(0.0, 0.0, 0.0, 1.0)]), 100)
    #coords = [(1:s).*25 for s ∈ size(vol)] # 25 μm resolution
    #vol = reverse(permutedims(vol, [1 3 2]), dims=(1,3))
    #coords = [1:x for x ∈ size(vol)].*50
    s = Scene(backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.0), resolution = resolution)
    coords = [1:s for s ∈ size(rotatereferenceatlas(vol))./2]
    coords[3] = .-coords[3]

    markercolor = :gray
    if ids isa Symbol
        if ids == :all
            # Get and plot a bunch of structures
            #ids = AN.getallstructureids()
            anns = unique(["FRP", "MO", "SS", "GU", "VISC", "AUD", "ACA", "PL", "ILA", "ORB", "AI", "RSP", "PTLp", "TEa", "PERI", "ECT", "OLF", "VISp", "VISl", "VISrl", "VISam", "VISpm", "VIS", "VISal","VISmma","VISmmp","VISli", "LGd","LD", "LP", "VPM", "TH", "MGm","MGv","MGd","PO","LGv","VL",
            "VPL","POL","Eth","PoT","PP","PIL","IntG","IGL","SGN","VPL","PF","RT", "CA1", "CA2","CA3", "DG", "SUB", "POST","PRE","ProS","HPF", "MB","SCig","SCiw","SCsg","SCzo","PPT","APN","NOT","MRN","OP","LT","RPF","CP"])
        elseif ids == :corticaltargets
            anns = unique(["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"])
        elseif ids == :targets # As listed in the visual coding white paper
            anns =  unique(["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam", "CA1", "CA3", "DG", "SUB", "ProS", "LGd", "LP", "APN"])
        else
            @error "`ids` keyword argument not valid "
        end
        t = AN.getstructuretree()
        D = t.get_id_acronym_map()
        ids = getindex.((D,), anns)
    end

    f = volume!(s, coords..., rotatereferenceatlas(vol)[1:2:end, 1:2:end, 1:2:end]; algorithm=:mip, colorrange=extrema(vol), colormap=collect(grad), ticks=nothing, fxaa=true, kwargs...)
    #s.plots[1].attributes[:fxaa] = true

    chrome = FileIO.load(download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/E6BF3C_5A4719_977726_FCFC82.png"))
    for probeid in unique(channels.probe_id)
        (x, y, z) = AN.getprobecoordinates(S, probeid)./50
        if !any(length.((x, y, z)).==0)
            f = meshscatter!(s, x, z, -y; markersize=1.0, fxaa=true, color=markercolor, matcap=chrome, shading=true, kwargs...)

            if dotext
                _, idx = findmax(z) # The highest unit
                text!(s, string(probeid), position=Vec3f0(x[idx], z[idx], -y[idx]), space=:data, textsize=5, align=(:right, :bottom), rotation=Quaternion((-0.3390201, 0.33899, 0.6205722, -0.620517)))
            end
        end
    end
    #scale!(s, 1, 1, -1)
    if dostructures
        channels = subset(AN.getchannels(),    :ecephys_probe_id           => ByRow(!ismissing),
                                                :ecephys_structure_id       => ByRow(!ismissing),
                                                :ecephys_structure_acronym  => ByRow(x->x!=="grey"))
        if isnothing(ids)
            ids = unique(channels[!, :ecephys_structure_id])
            #ids = ids[AN.getstructuretreedepth.(ids) .< 9]
        end
        for id ∈ ids
            try
                try
                    mask, _ = AN.getstructuremask(id)
                catch
                    mask = AN.buildstructuremask(id)
                end
                mask = Array{Float64}(rotatereferenceatlas(mask)[1:2:end, 1:2:end, 1:2:end])
                mc_algo = NaiveSurfaceNets(iso=0, insidepositive=true)
                m = GeometryBasics.Mesh(mask, mc_algo; origin=[min(coords[i]...) for i ∈ 1:length(coords)], widths=[abs(-(extrema(coords[i])...)) for i ∈ 1:length(coords)])
                c = AN.getstructurecolor(id)
                f = mesh!(s, m; color=RGBA(c.r, c.g, c.b, 0.41), fxaa=true, shading=true, kwargs...)
            catch y
                @warn y
            end
        end
    end
    return s
end

function exportreferencevolume(S, file::String="plot.html"; ids=:targets, kwargs...)
    open(file, "w") do io
            println(io, """
            <html>
                <head>
                </head>
                <body>
            """)
            show(io, MIME"text/html"(), Page(exportable=true, offline=true))
            fig = plotreferencevolume(S; dostructures=true, ids, kwargs...)
            rotate_cam!(fig, Vec3f0(0, 2.35, 0))
            show(io, MIME"text/html"(), fig)
            println(io, """
                </body>
            </html>
            """)
    end
end



function formattedreferencevolume(S, file::String="plot.html")
    exportreferencevolume(S, file;
                                    dostructures=true,
                                    ids=:targets,
                                    show_axis=false,
                                    shading=true, resolution=(1080, 1080));
end
