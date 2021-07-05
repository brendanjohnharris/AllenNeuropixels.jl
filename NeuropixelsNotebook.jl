### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 98c9bbd2-aac5-4c90-ac0c-d8d935f5cdaf
# For WGLMakie
begin
	using Pkg
	Pkg.add("JSServe")
	using JSServe
	Page()
end

# ╔═╡ 0a622047-c238-49c3-bdf0-3248d0f0d261
# Just for me
begin
	using Revise
	cd("../")
	Pkg.activate("@v1.6")
	Pkg.develop(path="./AllenNeuropixels.jl")
	cd("./AllenNeuropixels.jl")
	Pkg.activate("./")
end

# ╔═╡ 2e3d037c-932d-4d5a-afdb-351e836bdfb2
using WGLMakie

# ╔═╡ 79b92c17-cc5f-4ca2-8b08-f6015729a9a9
using DataFrames

# ╔═╡ 1a3eaeb7-db19-4bc3-b61a-84dd9caad009
using Statistics

# ╔═╡ c1397e8e-3ac2-4f49-8d88-cc1400a3b93e
using FileIO

# ╔═╡ 766a8af5-4c89-4fe7-883d-d960ef91edfd
md"""
# Allen Neuropixels
_Investigating low-dimensional projections of the Allen neuropixels visual coding of dataset_
"""

# ╔═╡ 66677c85-5aed-4793-aa03-ab070aa42dd0
html"""<img src="https://dl.dropbox.com/s/8amsi9wv0xtc1p3/muwfg7yln8a61.jpg?dl=0" alt="Bloons">"""

# ╔═╡ bea5de79-6e8a-42d8-ab76-bae8e3c23747
md"""
## Introduction
........intro..........

Details on neuropixels and the visual coding dataset can be found in the Allen SDK [docs](https://allensdk.readthedocs.io/en/latest/visual_coding_neuropixels.html), the [white-paper](https://dl.dropbox.com/s/tav6rft6iyd173k/neuropixels_visual_coding_-_white_paper_v10.pdf) or the [cheat-sheet](https://dl.dropbox.com/s/9v1x9eibmtb8fwg/neuropixels_cheat_sheet_nov_2019.pdf)
"""

# ╔═╡ f9ac9f6e-f129-4542-81a8-36e6cef9857a
md"## Required packages"

# ╔═╡ c48bd78e-aab0-49c0-b137-567c208b8aa1
md"""
Interfacing with the Python [Allen SDK](https://github.com/AllenInstitute/AllenSDK) is handled by the [AllenNeuropixels](https://github.com/brendanjohnharris/AllenNeuropixels.jl) package, which uses the [Makie](https://github.com/JuliaPlots/Makie.jl) package for plotting (since, unlike [Plots](https://github.com/JuliaPlots/Plots.jl), it is geared towards interactivity and can leverage the GPU for big datasets). We will use the WGLMakie as the backend, as well as DataFrames and Statistics for working with tables.
"""

# ╔═╡ c445ccf4-cf10-43b9-9c01-4051abc400ba
# Pkg.add(url="https://github.com/brendanjohnharris/AllenNeuropixels.jl#main")

# ╔═╡ d8fd4541-08a5-4ace-a998-769771c976e8
import AllenNeuropixels as AN

# ╔═╡ 5bfaefae-11a5-4567-8b83-d583f03a75a8
md"""
## Choosing a session
The Allen neuropixels visual coding dataset is subdivided into sessions. A session contains data from one murine subject across a recording interval during which it was shown a variety of visual stimuli (such as natural scenes, drift gratings and gabors). This data includes local field potentials (LFP) and spike times from each of the 374 data channels on, usually, six neuropixels probes inserted around the visual cortex. The LFP data is also downsampled by a half in time and a quarter over channels.

The entire neuropixels visual coding dataset contains dozens of sessions and is hundreds of megabytes in size, so we will first choose one session (a few gigabytes of data) instead of performing group-level analyses. To produce a dataframe of session details:
"""

# ╔═╡ 754f94e0-ccb2-4dc0-a534-ae94dd02bc02
session_table = AN.getsessiontable()

# ╔═╡ a2f81f74-36ed-42d4-89d5-828327c67318
md"""We'll pick a session that has six probes:"""

# ╔═╡ efd3ef52-4acb-4ffd-b794-1dfc4d9819c8
function numprobes(sessionid)
	session_table[session_table[!, :id].==sessionid, :probe_count][1]
end

# ╔═╡ e30e349c-6ad7-402b-90d1-7720b85a9c2c
md"And look for sessions that have probes intersecting the target regions listed in the white paper:"

# ╔═╡ d2d20884-0fcd-4ac7-8a14-dbf936445c3b
function targetintersections(sessionid)
	s = session_table[session_table.id.==sessionid, :ecephys_structure_acronyms][1]
	targets = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam", "CA1", "CA3", "DG", "SUB", "ProS", "LGd", "LP", "APN"];
	structures = eachmatch(r"'(.*?)'", s)
	structures = [replace(x.match, r"'"=>s"") for x ∈ structures]
	targetsintersected = [x ∈ structures for x ∈ targets]
	return mean(targetsintersected)
end;

# ╔═╡ 99fe41d6-661d-4b9f-a853-2f32ace53d72
md"""
We can filter down sessions in a few more ways. Firstly, the visaul coding data is divided into two stimulus sets. To summuraise the white paper:
- The **Brain Observatory 1.1** stimulus set:
    - Gabor patches appearing randomly on a 9 x 9 grid 
    - Dark or light flashes
    - Drift gratings in 8 directions and 5 temporal frequencies, with 15 repeats per condition
    - Static gratings at 6 orientations, 5 spatial frequencies, and 4 phases 
    - 118 natural images
    - Two natural movies from Touch of Evil; a 30 second clips repeated 20 times and a 120 second clip repeated 10 times


- The **Functional Connectivity** stimulus set:
    - Gabor patches
    - Dark or light flashes
    - Drift gratings in 4 directions and one temporal frequency with 75-80 repeats
    - Drift gratings in 4 directions and 9 contrasts
    - One natural movie shown 60 times plus 20 repeats of a temporally shuffled version
    - A dot motion stimulus of approximately 200 white dots on a gray background moving at one of 7 speeds in four different directions

In summary, the Brain Observatory dataset is formed by a greater variety of stimuli but a smaller number of repeats than the Functional Connectivity dataset. **We choose the Brain Observatory dataset since differences in activity patterns are likely to be more pronounced between different types of stimuli than variations on a stimulus, and our uniformly windowed analysis is not concerned with averaging over repetitions.**
"""

# ╔═╡ 82ce6d0f-b60b-41f2-bdce-c6ecf545bf65
md"""
Next, we can inspect the unit quality metrics of each session. Three metric criteria have recommended values:
- `amplitude_cutoff_maximum = 0.1`: Limits the approximate false negative rate in spike detection, calculated by assuming the distribution of amplitudes is symmetric. 
- `presence_ratio_minimum = 0.9`: Ensures proportion of 100 recording sub-intervals that have at least one spike detection is above 90% (i.e. the unit has not drifted and the sorting algorithm is correctly identifying spikes). 
- `isi_violations_maximum = 0.5`: Limits the number of units that record from multiple neurons. Inter-spike interval (ISI) violations are detections of spikes that occur during the typical refractory period of a single neuron, so are likely to originate from a second neuron.


To access a dataframe of metrics for the brain observatory dataset with the default metric filters:
"""

# ╔═╡ 5d47cb91-d8c7-41df-9778-c9a77266ba93
# This can take a while
metrics = AN.getunitanalysismetricsbysessiontype("brain_observatory_1.1")

# ╔═╡ 7d00cd72-a2c0-45a9-b777-b347286f7390
md"Some session-level metrics are:"

# ╔═╡ fd03139c-4bd1-4f26-a8c3-46c3feefd9c5
session_metrics = combine(groupby(metrics, :ecephys_session_id),
		:ecephys_session_id => numprobes∘unique => :num_probes,
		:ecephys_session_id => targetintersections∘unique => :target_intersections,
		:has_lfp_data=>all,
		:genotype => (x->all(isequal.(x, ("wt/wt",)))) => :is_wt,
		:max_drift=>median∘skipmissing,
		:d_prime=>median∘skipmissing,
		:isolation_distance=>median∘skipmissing,
		:silhouette_score=>median∘skipmissing,
		:snr=>median∘skipmissing)

# ╔═╡ 51eccd32-d4e9-4330-99f5-fca0f8534e43
md""""
Of particular importance are the first five columns. The ideal session will have 6probes, LFP data for all units and for now we wish to choose a wild type mouse. It should also contain data for probes that intersect with the major regions of interest highlighted in the white paper. The maximum drift, which measures mean (and possibly variance) nonstationarity unrelated to the visual stimulus (such as probe motion relative to the brain) should also be minimised. We can also look for probes that intersect all of the target regions, given in the white paper as:
"""

# ╔═╡ 927605f4-0b59-4871-a13f-420aadedd487
oursession = subset(session_metrics, 
						:num_probes => ByRow(==(6)),
						:target_intersections => ByRow(>(0.9)),
						:is_wt => ByRow(==(true)),
						:has_lfp_data_all => ByRow(==(true)),
						:max_drift_median_skipmissing => ByRow(<(25.0)))

# ╔═╡ 86382404-f784-4357-ae3c-05f7246d093b
md"""
To create a variable of type Session, which emulates Python's OOP interface:
"""

# ╔═╡ c3093ce3-7b73-49d4-8ce8-aaea4b49b685
session = AN.Session(oursession.ecephys_session_id[1])

# ╔═╡ 12d66773-c567-4e5a-a119-4fa5e06ec98c
md"We are only interested in channels located within the brain, so the channels for this session are:"

# ╔═╡ 5f944e90-5232-4508-9a3a-c678e6804104
channels = subset(AN.getchannels(session), :ecephys_structure_id=>ByRow(!ismissing))

# ╔═╡ 423edd0f-60ed-4d9a-b84a-47e47e560ae2
md"""
We can then plot the probe locations on the reference atlas. Note that the colored regions correspond to the target structures, and the straight probes have been warped into the common coordinate framework:
"""

# ╔═╡ 12f8e03b-4ea3-4211-a9b1-f8432dfae3a9
#s = AN.Plots.plotreferencevolume(session; dostructures=true,
#								ids=:targets,
#								show_axis=false,
#								shading=true); rotate_cam!(s, Vec3f0(0, 2.35, 0)); s

# To save you from waiting for the structure masks to download:
@html_str read(download("https://dl.dropbox.com/s/se2doygr56ox8hs/probelocations.html?dl=0"), String) 
# This renders best in firefox

# ╔═╡ Cell order:
# ╠═98c9bbd2-aac5-4c90-ac0c-d8d935f5cdaf
# ╠═0a622047-c238-49c3-bdf0-3248d0f0d261
# ╟─766a8af5-4c89-4fe7-883d-d960ef91edfd
# ╟─66677c85-5aed-4793-aa03-ab070aa42dd0
# ╟─bea5de79-6e8a-42d8-ab76-bae8e3c23747
# ╟─f9ac9f6e-f129-4542-81a8-36e6cef9857a
# ╟─c48bd78e-aab0-49c0-b137-567c208b8aa1
# ╠═c445ccf4-cf10-43b9-9c01-4051abc400ba
# ╠═d8fd4541-08a5-4ace-a998-769771c976e8
# ╠═2e3d037c-932d-4d5a-afdb-351e836bdfb2
# ╠═79b92c17-cc5f-4ca2-8b08-f6015729a9a9
# ╠═1a3eaeb7-db19-4bc3-b61a-84dd9caad009
# ╠═c1397e8e-3ac2-4f49-8d88-cc1400a3b93e
# ╟─5bfaefae-11a5-4567-8b83-d583f03a75a8
# ╠═754f94e0-ccb2-4dc0-a534-ae94dd02bc02
# ╟─a2f81f74-36ed-42d4-89d5-828327c67318
# ╠═efd3ef52-4acb-4ffd-b794-1dfc4d9819c8
# ╟─e30e349c-6ad7-402b-90d1-7720b85a9c2c
# ╠═d2d20884-0fcd-4ac7-8a14-dbf936445c3b
# ╟─99fe41d6-661d-4b9f-a853-2f32ace53d72
# ╟─82ce6d0f-b60b-41f2-bdce-c6ecf545bf65
# ╠═5d47cb91-d8c7-41df-9778-c9a77266ba93
# ╟─7d00cd72-a2c0-45a9-b777-b347286f7390
# ╠═fd03139c-4bd1-4f26-a8c3-46c3feefd9c5
# ╟─51eccd32-d4e9-4330-99f5-fca0f8534e43
# ╠═927605f4-0b59-4871-a13f-420aadedd487
# ╟─86382404-f784-4357-ae3c-05f7246d093b
# ╠═c3093ce3-7b73-49d4-8ce8-aaea4b49b685
# ╟─12d66773-c567-4e5a-a119-4fa5e06ec98c
# ╠═5f944e90-5232-4508-9a3a-c678e6804104
# ╟─423edd0f-60ed-4d9a-b84a-47e47e560ae2
# ╠═12f8e03b-4ea3-4211-a9b1-f8432dfae3a9
