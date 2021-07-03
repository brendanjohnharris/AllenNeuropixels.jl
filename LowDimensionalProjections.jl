### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 98c9bbd2-aac5-4c90-ac0c-d8d935f5cdaf
begin
	using Pkg
	Pkg.add("JSServe")
	using JSServe
	Page()
end

# ╔═╡ c445ccf4-cf10-43b9-9c01-4051abc400ba
Pkg.add(url="https://github.com/brendanjohnharris/AllenNeuropixels.jl#main")

# ╔═╡ 26dc1329-928e-4ca9-b8e7-beb3ddb05026
Pkg.add("WGLMakie")

# ╔═╡ de92ba74-dac7-4788-80f2-55f9f525c5d9
Pkg.add("PlutoUI")

# ╔═╡ 2e3d037c-932d-4d5a-afdb-351e836bdfb2
using WGLMakie

# ╔═╡ 2f2b7f62-a8bc-44f6-a9fa-aa13531d11e8
using PlutoUI

# ╔═╡ 766a8af5-4c89-4fe7-883d-d960ef91edfd
md"""
# Allen Neuropixels
_Investigating low dimensional projections of the Allen neuropixels visual coding of dataset_
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
Interfacing with the Python [Allen SDK](https://github.com/AllenInstitute/AllenSDK) is handled by the [AllenNeuropixels](https://github.com/brendanjohnharris/AllenNeuropixels.jl) package, which uses the [Makie](https://github.com/JuliaPlots/Makie.jl) package for plotting (since, unlike [Plots](https://github.com/JuliaPlots/Plots.jl), it is geared towards interactivity and can leverage the GPU for big datasets). We will use the WGLMakie as the backend.
"""

# ╔═╡ d8fd4541-08a5-4ace-a998-769771c976e8
import AllenNeuropixels as AN

# ╔═╡ 5bfaefae-11a5-4567-8b83-d583f03a75a8
md"""
## Choosing a session
The Allen neuropixels visual coding dataset is subdivided into sessions. A session contains data from one murine subject across a recording interval during which it was shown a variety of visual stimuli (such as natural scenes, drift gratings and gabors). This data includes local field potentials (LFP) and spike times from each of the 374 channels on, usually, six neuropixels probes inserted around the visual cortex. The LFP data is also downsampled by a half in time and a quarter over channels.

The entire neuropixels visual coding dataset contains dozens of sessions and is hundreds of megabytes in size, so we will first choose one session (a few gigabytes of data) instead of performing group-level analyses. To produce a dataframe of session details:
"""

# ╔═╡ 754f94e0-ccb2-4dc0-a534-ae94dd02bc02
𝒮 = AN.getsessiontable()

# ╔═╡ ba7ba800-0089-4121-98c8-4a8e3844a9c4
PlutoUI.WithIOContext(𝒮, displaysize=(9999,9999), :compact => true)

# ╔═╡ 94fb27ab-82fe-4f9e-a888-92f75d370f66
html"""
<iframe src="https://gistcdn.githack.com/brendanjohnharris/6e3c9a534898705ffdeeb8767d9fc25e/raw/e082f4c5deeb94b22e800369a021684b2980d9b6/plot.html"></iframe>
"""

# ╔═╡ e5be98e9-c57e-48a6-bfa1-8c871f2efd51
x = randn(100); y = x.^3; scatter(x, y)

# ╔═╡ 6098a444-4e8d-41e3-8fb1-1b6ba76afe72
lines(x, y)

# ╔═╡ a9fbd8e0-b1a5-4a2d-bdb9-b6909891e50f
pwd()

# ╔═╡ Cell order:
# ╠═98c9bbd2-aac5-4c90-ac0c-d8d935f5cdaf
# ╟─766a8af5-4c89-4fe7-883d-d960ef91edfd
# ╠═66677c85-5aed-4793-aa03-ab070aa42dd0
# ╠═bea5de79-6e8a-42d8-ab76-bae8e3c23747
# ╟─f9ac9f6e-f129-4542-81a8-36e6cef9857a
# ╠═c48bd78e-aab0-49c0-b137-567c208b8aa1
# ╠═c445ccf4-cf10-43b9-9c01-4051abc400ba
# ╠═26dc1329-928e-4ca9-b8e7-beb3ddb05026
# ╠═d8fd4541-08a5-4ace-a998-769771c976e8
# ╠═2e3d037c-932d-4d5a-afdb-351e836bdfb2
# ╠═5bfaefae-11a5-4567-8b83-d583f03a75a8
# ╠═754f94e0-ccb2-4dc0-a534-ae94dd02bc02
# ╠═de92ba74-dac7-4788-80f2-55f9f525c5d9
# ╠═2f2b7f62-a8bc-44f6-a9fa-aa13531d11e8
# ╠═ba7ba800-0089-4121-98c8-4a8e3844a9c4
# ╠═94fb27ab-82fe-4f9e-a888-92f75d370f66
# ╠═e5be98e9-c57e-48a6-bfa1-8c871f2efd51
# ╠═6098a444-4e8d-41e3-8fb1-1b6ba76afe72
# ╠═a9fbd8e0-b1a5-4a2d-bdb9-b6909891e50f
