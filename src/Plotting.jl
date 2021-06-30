using WGLMakie
using JSServe
using Colors
using CartesianBoxes

# This works best from Julia in the terminal, not in VSCode

Page(exportable=true, offline=true)
rotatereferenceatlas = x -> reverse(permutedims(x, (1, 3, 2)), dims=(3,))
# Plot the reference volume and probe locations
function plotreferencevolume(S; dostructures = true, resolution=(1200, 800))
    channels = getchannels(S)
    vol, info = gettemplatevolume()
    vol = Array{Float16}(vol)
    grad = cgrad(cgrad([RGBA(1.0, 1.0, 1.0, 0.0), RGBA(0.0, 0.0, 0.0, 1.0)]), 100)
    #coords = [(1:s).*25 for s ∈ size(vol)] # 25 μm resolution
    #vol = reverse(permutedims(vol, [1 3 2]), dims=(1,3))
    #coords = [1:x for x ∈ size(vol)].*50
    s = Scene()
    coords = [1:s for s ∈ size(rotatereferenceatlas(vol))./2]
    coords[3] = .-coords[3]
    f = volume!(s, coords..., rotatereferenceatlas(vol)[1:2:end, 1:2:end, 1:2:end]; algorithm=:mip, colorrange=extrema(vol), colormap=collect(grad), figure =(resolution = resolution,), ticks=nothing) # Or :mip

    if dostructures
        channels = subset(get_channels(),       :ecephys_probe_id           => ByRow(!ismissing),
                                                :ecephys_structure_id       => ByRow(!ismissing),
                                                :ecephys_structure_acronym  => ByRow(x->x!=="grey"))

        ids = unique(channels[!, :ecephys_structure_id])
        println(getstructuretreedepth.(ids))
        ids = ids[getstructuretreedepth.(ids) .< 9]
        for id ∈ ids[5:10]
            try
                mask, _ = getstructuremask(id)
                mask = Array{Float16}(mask)
                idxs = CartesianBoxes.boundingbox(rotatereferenceatlas(mask)[1:2:end, 1:2:end, 1:2:end])
                tidxs = Tuple.(CartesianIndices(idxs))
                xx, yy, zz = [unique([i[j] for i ∈ tidxs]) for j ∈ 1:3]
                c = getstructurecolor(id)
                maskgrad = :transparent => RGBA(c.r, c.b, c.g, 0.5)
                volume!(s, coords[1][xx], coords[2][yy], coords[3][zz], rotatereferenceatlas(mask)[1:2:end, 1:2:end, 1:2:end][idxs]; algorithm=:mip, colorrange=extrema(mask), color=(c, 0.5), transparency=true)
            catch y
                @warn y
            end
        end

    end

    (x, y, z) = getprobecoordinates(S)./50
    meshscatter!(s, x, z, -y, markersize=1.0)
    #scale!(s, 1, 1, -1)
    return s
end

function exportreferencevolume(S, file::String="plot.html")
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
        show(io, MIME"text/html"(), plotreferencevolume(S, resolution=(2400, 1600)))
        println(io, """
            </body>
        </html>
        """)
    end
end