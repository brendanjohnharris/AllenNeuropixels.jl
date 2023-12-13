using Makie
using WGLMakie
using FileIO
using Makie.Colors
import AllenNeuropixels as AN
import AllenNeuropixels.AllenNeuropixelsBase as ANB

S = ANB.VisualBehavior.Session(1044385384)
dotext = true
dostructures = true
ids = :targets
size = (1920, 1080)

channels = AN.getchannels(S)
vol, info = ANB.gettemplatevolume()
vol = Array{Float16}(vol)
grad = cgrad(cgrad([RGBA(1.0, 1.0, 1.0, 0.0), RGBA(0.0, 0.0, 0.0, 1.0)]), 100)
#coords = [(1:s).*25 for s ∈ size(vol)] # 25 μm size
#vol = reverse(permutedims(vol, [1 3 2]), dims=(1,3))
#coords = [1:x for x ∈ size(vol)].*50
# s = Scene(backgroundcolor = RGBA(1.0, 1.0, 1.0, 0.0), size = size, transparency=true)
coords = [1:s for s in size(AN.Plots.rotatereferenceatlas(vol)) ./ 2]
coords[3] = .-coords[3]

## RPR stuff
radiance = 500
lights = [EnvironmentLight(1.0, load(RPR.assetpath("studio032.exr"))),
    PointLight(Vec3f(10), RGBf(radiance, radiance, radiance * 1.1))]

fig = Figure();
radiance = 10000
lights = [
    EnvironmentLight(1, Makie.FileIO.load(RPR.assetpath("studio027.exr"))),
    PointLight(Vec3f(10, 10, 10), RGBf(radiance, radiance, radiance)),
]
ax = LScene(fig[1, 1]; show_axis = false, scenekw = (lights = lights,))
screen = RPRMakie.Screen(ax.scene; iterations = 10, plugin = RPR.Northstar);
matsys = screen.matsys;
context = screen.context;
mat = RPR.Plastic(matsys)
cam = cameracontrols(ax.scene)
cam.eyeposition[] = Vec3f(1.5, 0, 0)
cam.lookat[] = Vec3f(0, 0, 0)
cam.upvector[] = Vec3f(0, 0, 1)
cam.fov[] = 50
mesh!(ax, Sphere(Point3f(0), 1), material = mat)
screen

# image = colorbuffer(screen)
# save("materials.png", image)

markercolor = :gray
if ids isa Symbol
    if ids == :all
        # Get and plot a bunch of structures
        #ids = AN.getallstructureids()
        anns = unique(["FRP", "MO", "SS", "GU", "VISC", "AUD", "ACA", "PL", "ILA", "ORB",
                          "AI", "RSP", "PTLp", "TEa", "PERI", "ECT", "OLF", "VISp", "VISl",
                          "VISrl", "VISam", "VISpm", "VIS", "VISal", "VISmma", "VISmmp",
                          "VISli", "LGd", "LD", "LP", "VPM", "TH", "MGm", "MGv", "MGd",
                          "PO", "LGv", "VL",
                          "VPL", "POL", "Eth", "PoT", "PP", "PIL", "IntG", "IGL", "SGN",
                          "VPL", "PF", "RT", "CA1", "CA2", "CA3", "DG", "SUB", "POST",
                          "PRE", "ProS", "HPF", "MB", "SCig", "SCiw", "SCsg", "SCzo", "PPT",
                          "APN", "NOT", "MRN", "OP", "LT", "RPF", "CP"])
    elseif ids == :corticaltargets
        anns = unique(["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam"])
    elseif ids == :targets # As listed in the visual coding white paper
        anns = unique([
                          "VISp",
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
                          "APN",
                      ])
    else
        @error "`ids` keyword argument not valid "
    end
    t = ANB.getstructuretree() # Fucked
    D = t.get_id_acronym_map() |> Dict
    ids = getindex.((D,), anns)
end

f = volume!(ax, coords..., AN.Plots.rotatereferenceatlas(vol)[1:2:end, 1:2:end, 1:2:end];
            algorithm = :mip, colorrange = extrema(vol), colormap = collect(grad),
            ticks = nothing, fxaa = true, transparency = true)
#s.plots[1].attributes[:fxaa] = true

chrome = FileIO.load(download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/E6BF3C_5A4719_977726_FCFC82.png"))
for probeid in unique(channels.probe_id)
    (x, y, z) = AN.getprobecoordinates(S, probeid) ./ 50
    if !any(length.((x, y, z)) .== 0)
        f = meshscatter!(s, x, z, -y; markersize = 1.0, fxaa = true, color = markercolor,
                         matcap = chrome, shading = true, transparency = true, kwargs...)

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
            mc_algo = NaiveSurfaceNets(iso = 0, insidepositive = true)
            m = GeometryBasics.Mesh(mask, mc_algo;
                                    origin = [min(coords[i]...) for i in 1:length(coords)],
                                    widths = [abs(-(extrema(coords[i])...))
                                              for i in 1:length(coords)])
            c = AN.getstructurecolor(id)
            f = mesh!(s, m; color = RGBA(c.r, c.g, c.b, 0.41), fxaa = true, shading = true,
                      transparency = true, kwargs...)
        catch y
            @warn y
        end
    end
end
