using GLMakie
GLMakie.activate!(ssao = true)
# using CairoMakie
using FileIO
import AllenNeuropixels as AN
import AllenNeuropixels.AllenNeuropixelsBase as ANB
using DataFrames
import AllenNeuropixelsBase.GeometryBasics
using Normalization
using DimensionalData
import AllenNeuropixels.Plots.plotbrain!
using Foresight
set_theme!(foresight())

LFP = load(abspath(Base.find_package("AllenNeuropixels"), "..", "..", "test",
                   "PlotBrain.jld2"))["LFP"]
sessionid = LFP[1].metadata[:sessionid]
channels = dims.(LFP, 2) .|> collect

S = ANB.VisualBehavior.Session(sessionid)

f = Figure(; resolution = (1080, 720), lightposition = Vec3(100000.0, 100000.0, 0.0))
ax = Axis3(f[1, 1]; aspect = :data)
hidedecorations!(ax)
ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
probeids, c, p = AN.Plots.plotbrain!(ax, S; dark = false, probestyle = :meshscatter,
                                     channels, markersize = 100.0, fontsize = 15.0)
f

# ? Animate activity from an LFP trace. An example is saved in PlotBrain.jld2
framerate = 60
ts = size(LFP[1], 1)
(idxs = indexin(probeids, getindex.(DimensionalData.metadata.(LFP), :probeid)); X = collect.(LFP[idxs])) # Sort to probe order in plot

Nc = length.([_c[] for _c in c])

colormap = cgrad(:binary)
X = [MinMax(x, dims = 1)(x) for x in X]
Nt = unique(size.(X, 1))
length(Nt) > 1 && error("All provided arrays must have the same number of rows")
Nt = first(Nt)
record(f, "PlotBrain.gif", 1:4:Nt; framerate) do t
    @info "Rendering $t of $Nt frames"
    for i in eachindex(c)
        c[i][] = colormap[X[i][t, :]]
    end
    if t < Nt ÷ 2
        ax.azimuth[] += 0.001
        ax.elevation[] += 0.0005
    else
        ax.azimuth[] -= 0.001
        ax.elevation[] -= 0.0005
    end
end

# ? Do it in a loop for GLMakie interactivity
if false
    map(ts) do t
        for i in eachindex(c)
            c[i][] = colormap[X[i][t, :]]
        end
    end
    sleep(1 / fs)
end
