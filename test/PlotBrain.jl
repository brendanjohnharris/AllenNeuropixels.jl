@info @__FILE__
@info pwd()
using IntervalSets
using Normalization
using GLMakie
import GLMakie.record
using DimensionalData
using ProgressLogging
using Foresight
using TimeseriesTools
using AllenNeuropixels
import AllenNeuropixels as AN
set_theme!(foresight(:dark))

plotpath = "/home/brendan/OneDrive/Masters/Code/Vortices/Julia/AllenAttention/test/Neuropixels/VisualBehavior/SingleSubject/ThetaBurstBehavior/"

# * Get data
begin
    sessionid = 1130113579 # 1067588044 #
    pass = (3, 10)
    fs = 1250
end
begin
    session = AN.Session(sessionid)
    probeids = AN.getprobes(session).id
    _LFP = map(x -> AN.getlfp(session, x; times = 4600 .. 4700), probeids)
end
begin
    start = 100
    step = 2
    N = 2000
    LFP = AN.rectifytime.(_LFP)
    LFP = bandpass.(LFP, [pass])
    LFP = getindex.(LFP, [range(; start, step, length = N)], [:])
    channels = dims.(LFP, 2) .|> collect
end

begin # Set up plot
    f = Figure(; size = (720, 1080), light_direction = Vec3(100000.0, 100000.0, 0.0))
    ax = Axis3(f[1, 1]; aspect = :data, title = "0 s", titlegap = -100)
    hidedecorations!(ax)
    ax.xspinesvisible = ax.yspinesvisible = ax.zspinesvisible = false
    probeids, c, p = AN.Plots.plotbrain!(ax, session; dark = false,
                                         probestyle = :meshscatter,
                                         channels, markersize = 100.0, fontsize = 15.0)
    f
end

begin # ? Animate activity from an LFP trace. An example is saved in PlotBrain.jld2
    framerate = 60
    ts = size(LFP[1], 1)
    truets = times(LFP[1])
    (idxs = indexin(probeids, getindex.(DimensionalData.metadata.(LFP), :probeid)); X = collect.(LFP[idxs])) # Sort to probe order in plot

    Nc = length.([_c[] for _c in c])

    colormap = cgrad(:binary)
    X = [MinMax(x, dims = 1)(x) for x in X]
    Nt = unique(size.(X, 1))
    length(Nt) > 1 && error("All provided arrays must have the same number of rows")
    Nt = first(Nt)
    @withprogress record(f, "PlotBrain.mp4", 1:Nt; framerate) do t
        for i in eachindex(c)
            c[i][] = colormap[X[i][t, :]]
        end
        if t < Nt ÷ 2
            ax.azimuth[] += 0.0015
            ax.elevation[] += 0.0005
        else
            ax.azimuth[] -= 0.0015
            ax.elevation[] -= 0.0005
        end
        ax.title[] = "$(round(truets[t], digits=1)) s"
        @logprogress t / Nt
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
