using CairoMakie
import AllenNeuropixels as AN
using DimensionalData

# * Choose a good session
sessionid = 757216464
session = AN.Session(sessionid)

# * Pick a probe
probeid = 769322749

# * Load LFP data for a single probe from VISp
fs = 1250
X = AN.getlfp(session, probeid; inbrain=200)

# * Plot psd for a single probe
fig1 = AN.Plots.plotLFPspectra(session, probeid, X)
current_axis().xscale=Makie.pseudolog10; current_axis().xticks = vcat(0:5:10, 20:20:100)
xlims!(current_axis(), [0, 100]); ylims!(current_axis(), [1e-4, 1])
# fig1 = AN.Plots.plotLFPspectra(session, probeid, X[1:10:end, :]);

# * Then partition out just the visp channels
channels = dims(X, :channel) |> collect
structures = AN.getstructureacronyms(channels)
_X = X[:, structures.=="VISp"]
fig2 = AN.Plots.plotLFPspectra(session, probeid, _X)
current_axis().xscale=Makie.pseudolog10; current_axis().xticks = vcat(0:5:10, 20:20:100)
xlims!(current_axis(), [0, 100]); ylims!(current_axis(), [1e-4, 1])

# * For low-contrast stimuli
ts = AN.getstimulustimes(session, "gabors")
times = [any(in.(t, ts)) for t in dims(X, Ti)] # Times in any of the stimulus epochs
_X = X[times, structures.=="VISp"]
fig3 = AN.Plots.plotLFPspectra(session, probeid, _X)
current_axis().xscale=Makie.pseudolog10; current_axis().xticks = vcat(0:5:10, 20:20:100)
xlims!(current_axis(), [0, 100]); ylims!(current_axis(), [1e-4, 1])

# * And finally for high-contrast stimuli
ts = AN.getstimulustimes(session, "static_gratings")
times = [any(in.(t, ts)) for t in dims(X, Ti)] # Times in any of the stimulus epochs
_X = X[times, structures.=="VISp"]
fig4 = AN.Plots.plotLFPspectra(session, probeid, _X)
current_axis().xscale=Makie.pseudolog10; current_axis().xticks = vcat(0:5:10, 20:20:100)
xlims!(current_axis(), [0, 100]); ylims!(current_axis(), [1e-4, 1])
