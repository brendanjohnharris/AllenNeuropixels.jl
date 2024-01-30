@info @__FILE__
@info pwd()
using GLMakie
using GLMakie.FileIO
using Downloads
using Foresight
using AllenNeuropixels
import AllenNeuropixels as AN
import AllenNeuropixelsBase as ANB
set_theme!(foresight())

begin
    sessionid = 1128520325
    session = AN.Session(sessionid)
    probeids = AN.getprobes(session).id
    chrome = FileIO.load(Downloads.download("https://raw.githubusercontent.com/nidorx/matcaps/master/1024/C8D1DC_575B62_818892_6E747B.png"))
end

begin # Set up plot
    f = Figure(; size = (2560, 1440))
    ax = Axis3(f[1, 1];
               aspect = :data,
               xzpanelvisible = false,
               yzpanelvisible = false,
               xypanelvisible = false)
    ax.scene.lights = [DirectionalLight(RGBf(8, 8, 8), Vec3f(2, -1, 1))]
    # ax.scene.lights[2] = PointLight(Makie.RGB(1.0, 1.0, 1),
    #                                       Vec3(14530.263, 11445.357, 9516.364))
    hidedecorations!(ax)
    ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
    structurecolors = cgrad(:inferno)[(1:6) ./ 7] .* 1.3
    # push!(structurecolors, Foresight.keppel)
    probeids, c, p = AN.Plots.plotbrain!(ax, session;
                                         dark = false,
                                         probestyle = :meshscatter,
                                         markersize = 100.0,
                                         fontsize = 15.0,
                                         matcap = false,
                                         shading = Makie.MultiLightShading)
    #  ids = [
    #      "VISp",
    #      "VISl",
    #      "VISrl",
    #      "VISal",
    #      "VISpm",
    #      "VISam",
    #  ]),
    #  structurecolors,
    #  meshparams = (; alpha = 0.8))
    # * Color each probe by its position in the cortical hierarchy...
    # meshscatter!(ax, [0], [0], [0], markersize = 2000.0, color = :black)
    ax.azimuth = 2.6
    ax.elevation = 0.24
    s = ["VISp", "VISl", "VISal", "VISrl", "VISpm", "VISam"]
    ss = AN.getprobestructures(session, s)
    ss = getindex.([ss], probeids)
    is = indexin(ss, s)

    # n = 100
    # cs = [vcat(fill(Makie.colorant"gold", length(c[i][]) - n),
    #            fill(structurecolors[is[i]], n)) for i in eachindex(is)]
    # [c[i][] = cs[i] for i in eachindex(is)]
    [c[i][] = fill(structurecolors[is[i]], length(c[i][])) for i in eachindex(is)]
    # Camera3D(ax.scene)
    f
    save("$(@__DIR__)/plotbrain_static.png", f)
end
