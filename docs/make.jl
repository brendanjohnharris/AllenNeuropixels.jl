using CairoMakie
using AllenNeuropixels
using Documenter

# DocMeta.setdocmeta!(AllenNeuropixels, :DocTestSetup, :(using CairoMakie; using AllenNeuropixels); recursive=true)

makedocs(;
         modules = [AllenNeuropixels],
         authors = "brendanjohnharris <brendanjohnharris@gmail.com> and contributors",
         repo = "https://github.com/brendanjohnharris/AllenNeuropixels.jl/blob/{commit}{path}#{line}",
         sitename = "AllenNeuropixels.jl",
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://brendanjohnharris.github.io/AllenNeuropixels.jl",
                                  edit_link = "main",
                                  assets = String[],),
         pages = [
             "Home" => "index.md",
         ],)

deploydocs(;
           repo = "github.com/brendanjohnharris/AllenNeuropixels.jl",
           devbranch = "main",)
