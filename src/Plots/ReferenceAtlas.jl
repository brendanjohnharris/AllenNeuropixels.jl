

rotatereferenceatlas = x -> reverse(permutedims(x, (1, 3, 2)), dims=(3,))
# Plot the reference volume and probe locations
function plotreferencevolume(S; dostructures = true, resolution=(2400, 1600), ids=nothing)
    #! Check probe locations are correct
    channels = AN.getchannels(S)
    vol, info = AN.gettemplatevolume()
    vol = Array{Float16}(vol)
    grad = cgrad(cgrad([RGBA(1.0, 1.0, 1.0, 0.0), RGBA(0.0, 0.0, 0.0, 1.0)]), 100)
    #coords = [(1:s).*25 for s ∈ size(vol)] # 25 μm resolution
    #vol = reverse(permutedims(vol, [1 3 2]), dims=(1,3))
    #coords = [1:x for x ∈ size(vol)].*50
    s = Scene()
    coords = [1:s for s ∈ size(rotatereferenceatlas(vol))./2]
    coords[3] = .-coords[3]
    f = volume!(s, coords..., rotatereferenceatlas(vol)[1:2:end, 1:2:end, 1:2:end]; algorithm=:mip, colorrange=extrema(vol), colormap=collect(grad), figure =(resolution = resolution,), ticks=nothing, fxaa=true) # Or :mip
    #s.plots[1].attributes[:fxaa] = true

    (x, y, z) = AN.getprobecoordinates(S)./50
    f = meshscatter!(s, x, z, -y, markersize=1.0, fxaa=true)

    #scale!(s, 1, 1, -1)
    if dostructures
        channels = subset(AN.get_channels(),       :ecephys_probe_id           => ByRow(!ismissing),
                                                :ecephys_structure_id       => ByRow(!ismissing),
                                                :ecephys_structure_acronym  => ByRow(x->x!=="grey"))
        if isnothing(ids)
            ids = unique(channels[!, :ecephys_structure_id])
            #ids = ids[AN.getstructuretreedepth.(ids) .< 9]
        end
        for id ∈ ids
            try
                mask, _ = AN.getstructuremask(id)
                mask = Array{Float64}(rotatereferenceatlas(mask)[1:2:end, 1:2:end, 1:2:end])
                mc_algo = NaiveSurfaceNets(iso=0, insidepositive=true)
                m = GeometryBasics.Mesh(mask, mc_algo; origin=[min(coords[i]...) for i ∈ 1:length(coords)], widths=[abs(-(extrema(coords[i])...)) for i ∈ 1:length(coords)])
                c = AN.getstructurecolor(id)
                f = mesh!(s, m; color=RGBA(c.r, c.b, c.g, 0.41), fxaa=true)
            catch y
                @warn y
            end
        end

    end

    return s
end

function exportreferencevolume(S, file::String="plot.html"; ids=nothing)
    open(file, "w") do io
        println(io, """
        <html>
            <head>
            </head>
            <body>
        """)
        # before doing anything else,
        # make sure the Page setup code gets rendered as HTML
        show(io, MIME"text/html"(), Page(exportable=true, offline=true))
        # Then, you can just inline plots or whatever you want :)
        show(io, MIME"text/html"(), plotreferencevolume(S; resolution=(2400, 1600), dostructures=true, ids))
        println(io, """
            </body>
        </html>
        """)
    end
end

