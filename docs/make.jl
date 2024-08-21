using Makie
using AllenNeuropixels
using Documenter

DocMeta.setdocmeta!(AllenNeuropixelsBase, :DocTestSetup, :(using AllenNeuropixelsBase);
                    recursive = true)
DocMeta.setdocmeta!(AllenNeuropixels, :DocTestSetup, :(using AllenNeuropixels);
                    recursive = true)

pages = ["Home" => "index.md",
    "AllenNeuropixels" => "allen_neuropixels.md",
    "AllenNeuropixelsBase" => "allen_neuropixels_base.md"]

makedocs(;
         modules = [AllenNeuropixels, AllenNeuropixelsBase],
         authors = "brendanjohnharris <bhar9988@uni.sydney.edu.au> and contributors",
         sitename = "AllenNeuropixels.jl",
         format = Documenter.HTML(;
                                  canonical = "https://brendanjohnharris.github.io/AllenNeuropixels.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages)

deploydocs(;
           repo = "github.com/brendanjohnharris/AllenNeuropixels.jl",
           devbranch = "main",)
