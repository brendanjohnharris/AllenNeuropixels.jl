# AllenNeuropixels

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://brendanjohnharris.github.io/AllenNeuropixels.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://brendanjohnharris.github.io/AllenNeuropixels.jl/dev/)
[![Build Status](https://github.com/brendanjohnharris/AllenNeuropixels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/brendanjohnharris/AllenNeuropixels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/brendanjohnharris/AllenNeuropixels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/brendanjohnharris/AllenNeuropixels.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13363714.svg)](https://doi.org/10.5281/zenodo.13363714)


This package contains high-level analysis and visualization tools for the Allen Neuropixels mouse electrophysiology datasets ([Visual Coding](https://portal.brain-map.org/circuits-behavior/visual-coding-neuropixels) and [Visual Behavior](https://portal.brain-map.org/circuits-behavior/visual-behavior-neuropixels)) in Julia.
It also re-exports low-level functions for accessing and formatting the data ([AllenNeuropixelsBase.jl](https://www.github.com/brendanjohnharris/AllenNeuropixelsBase.jl)) and bindings to the Python [AllenSDK](https://github.com/AllenInstitute/AllenSDK) ([AllenSDK.jl](https://www.github.com/brendanjohnharris/AllenSDK.jl)).


![PlotBrain](plotbrain.png)