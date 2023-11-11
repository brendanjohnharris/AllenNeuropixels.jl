# using GLMakie
# GLMakie.activate!(ssao = true)
using CairoMakie
using FileIO
import AllenNeuropixels as AN
import AllenNeuropixels.AllenNeuropixelsBase as ANB
using DataFrames
import AllenNeuropixelsBase.GeometryBasics

S = ANB.VisualBehavior.Session(1067588044)

f = Figure(; resolution = (1920, 1080))
ax = Axis3(f[1, 1]; aspect = :data)
hidedecorations!(ax)
ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
c, p = plotbrain!(ax, S; dark = false, probestyle = :meshscatter)
f

# ? Animate activity from an LFP trace. An example is saved in PlotBrain.jld2
LFP = joinpath(Base.find_package("AllenNeuropixels"), "test", "PlotBrain.jld2")
framerate = 30
timestamps = range(0, 2, step = 1 / framerate)

record(fig, "time_animation.mp4", timestamps;
       framerate = framerate) do t
    time[] = t
end
[_c[] = cgrad(:bone)[(1:length(_c[])) ./ length(_c[])] for _c in c]

# ? Do it in a loop for GLMakie interactivity
