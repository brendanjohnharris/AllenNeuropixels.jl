### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 6eecb9bc-167d-4730-8fb4-92454d0f42fc
# For WGLMakie
begin
	using Pkg
	Pkg.add("JSServe")
	using JSServe
	Page()
end

# ╔═╡ bd11a795-d634-4e6f-9df7-89ebd17dbaf8
# Just for me
begin
	using Revise
	cd("../")
	Pkg.activate("@v1.6")
	Pkg.develop(path="./AllenNeuropixels.jl")
	cd("./AllenNeuropixels.jl")
	Pkg.activate("./")
end

# ╔═╡ 63bceb60-9821-42f6-b71d-9838bd5a31d7
using WGLMakie

# ╔═╡ 89753cf2-e09a-47aa-a046-02d7cc7b7609
using DataFrames

# ╔═╡ 7ed3d5aa-f23b-4826-bae7-5293776dc44d
using Statistics

# ╔═╡ 97b11cb5-51c3-4e8b-a1fe-845a22c11161
using FileIO

# ╔═╡ 51fc0034-a12a-4734-a965-498cee9dc424
md"# Neuropixels Visual Coding Stimuli"

# ╔═╡ 5992f59a-d7d6-4dd7-9a12-19c1ed42d6e6
md"**Accessing and visualising the visual stimuli presented to mice for the Neuropixels visal coding sessions using the Allen SDK**"

# ╔═╡ 8c4a5164-e315-403e-8b9f-2985c8ad84ff
md"[yayaya](http://alleninstitute.github.io/AllenSDK/allensdk.brain_observatory.stimulus_info.html)"

# ╔═╡ Cell order:
# ╠═6eecb9bc-167d-4730-8fb4-92454d0f42fc
# ╠═bd11a795-d634-4e6f-9df7-89ebd17dbaf8
# ╠═63bceb60-9821-42f6-b71d-9838bd5a31d7
# ╠═89753cf2-e09a-47aa-a046-02d7cc7b7609
# ╠═7ed3d5aa-f23b-4826-bae7-5293776dc44d
# ╠═97b11cb5-51c3-4e8b-a1fe-845a22c11161
# ╠═51fc0034-a12a-4734-a965-498cee9dc424
# ╠═5992f59a-d7d6-4dd7-9a12-19c1ed42d6e6
# ╠═8c4a5164-e315-403e-8b9f-2985c8ad84ff
