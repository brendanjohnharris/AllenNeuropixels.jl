import AllenNeuropixels as AN
using DataFrames
using Test
using Random, Statistics
using Logging

@testset "AllenNeuropixels" begin
    cs = AN.getprobes()
    @test cs isa DataFrame
    @test !isempty(cs)
end

@testset "bandpass keeps the passband" begin
    fs = 1250.0
    ts = 0:(1 / fs):2
    x = AN.ToolsArray(sin.(2π * 50 .* ts) .+ sin.(2π * 5 .* ts) .+ sin.(2π * 300 .* ts), (AN.𝑡(ts),))
    y = AN.bandpass(x; pass = [30, 80])
    @test AN.dims(y) == AN.dims(x)
    mid = 500:2000 # away from the filter's edge transients
    @test parent(y)[mid] ≈ sin.(2π * 50 .* ts[mid]) atol = 0.05
end


# A synthetic LFP: drifting noise with 60 Hz bursts at 2, 5 and 8 s
const BURSTTIMES = (2.0, 5.0, 8.0)
function syntheticlfp(; fs = 1250.0, rng = MersenneTwister(1))
    ts = 0:(1 / fs):(10 - 1 / fs)
    x = cumsum(randn(rng, length(ts))) .* 0.02 .+ randn(rng, length(ts)) .* 0.2
    for t0 in BURSTTIMES
        x .+= 3.0 .* exp.(-((ts .- t0) .^ 2) ./ (2 * 0.04^2)) .* sin.(2π * 60 .* ts)
    end
    AN.ToolsArray(x, (AN.𝑡(ts),); refdims = (AN.Chan(1),))
end

@testset "detectbursts finds injected gamma bursts" begin
    res = AN.wavelettransform(syntheticlfp(); pass = (10, 200))
    @test res isa AN.WaveletMatrix
    B = AN.detectbursts(convert(AN.LogWaveletMatrix, res); pass = (30, 100))
    for t0 in BURSTTIMES
        b = B[argmin(abs.(AN.peaktime.(B) .- t0))]
        @test abs(AN.peaktime(b) - t0) < 0.1
        @test 50 < AN.peakfreq(b) < 70
    end

    ϕ = convert(AN.LogWaveletMatrix, angle.(AN._wavelettransform(syntheticlfp(); pass = (10, 200))))
    Bϕ = AN.addphasemasks(B, ϕ)
    @test all(b -> AN.phasemask(b) isa AN.LogWaveletMatrix, Bϕ)
    @test all(isnothing ∘ AN.phasemask, B) # the input bursts are left untouched
    m = AN.burstmask(convert(AN.LogWaveletMatrix, res), B[1].peak)
    @test m isa AN.LogWaveletMatrix
    @test size(m, AN.Log𝑓) < size(res, 2) # the mask is cut to the burst's band
end

@testset "every supra-threshold region becomes a burst" begin
    res = convert(AN.LogWaveletMatrix, AN.wavelettransform(syntheticlfp(); pass = (10, 200)))
    mask = (AN.burstthreshold(res, 4) .> 0) .& (AN.burstcurvature(res, 2) .> 0)
    boxes = AN.ImageMorphology.component_boxes(AN.ImageMorphology.label_components(Array(mask)))
    n = count(b -> all(size(b) .> 1), boxes[begin+1:end]) # regions one pixel thick are dropped by design
    @test n > 0
    @test length(AN._detectbursts(res)) == n
end

@testset "globalbursts keeps the first interval" begin
    res = convert(AN.LogWaveletMatrix, AN.wavelettransform(syntheticlfp(); pass = (10, 200)))
    B = AN.detectbursts(res; pass = (30, 100))
    B = filter(b -> any(t0 -> abs(AN.peaktime(b) - t0) < 0.1, BURSTTIMES), B)
    G = AN.globalbursts(Dict(1 => B, 2 => B))
    @test all(t0 -> any(g -> t0 ∈ g, G), BURSTTIMES)
end

@testset "fooofedwavelet" begin
    res = convert(AN.LogWaveletMatrix, AN.wavelettransform(syntheticlfp(); pass = (10, 200)))
    r = AN.fooofedwavelet(res)
    @test size(r) == size(res)
    @test all(isfinite, r)
end

@testset "wavelettransform keeps channel order" begin
    x = syntheticlfp()
    X = AN.ToolsArray(hcat(parent(x), 2 .* parent(x), 3 .* parent(x)),
                      (AN.dims(x, AN.𝑡), AN.Chan([3, 1, 2])))
    W = AN.wavelettransform(X; pass = (10, 200))
    w = AN.wavelettransform(x; pass = (10, 200))
    @test AN.lookup(W, AN.Chan) == [3, 1, 2]
    for (k, c) in enumerate([3, 1, 2]) # the magnitude scales with the amplitude
        @test parent(W[AN.Chan(AN.At(c))]) ≈ k .* parent(w)
    end
end

@testset "fitdiagonalgaussian" begin
    x = range(-1, 1, length = 41)
    y = range(0, 2, length = 31)
    p = [2.0, 0.2, 1.1, 0.3, 0.4]
    Z = [AN.gaussian2((a, b), p) for a in x, b in y]
    fit = @test_logs min_level = Logging.Warn AN.fitdiagonalgaussian(x, y, Z) # no deprecation warnings
    @test abs.(fit.param) ≈ p rtol = 1e-4
end

using CairoMakie # loads AllenNeuropixels.Plots

@testset "plots" begin
    res = convert(AN.LogWaveletMatrix, AN.wavelettransform(syntheticlfp(); pass = (10, 200)))
    @test AN.Plots.plotwaveletvariability(res, f -> 0.0) isa Figure

    B = AN.detectbursts(res; pass = (30, 100))
    @test AN.Plots.plotfit(B[1]) isa Union{Figure, Makie.FigureAxisPlot}

    # Neuropixels channels come in pairs at each depth, so depths are not a range
    layers = ["L1", "L1", "L2/3", "L2/3", "L2/3", "L2/3", "L4", "L4"]
    depths = [-100.0, -100.0, -80.0, -80.0, -60.0, -60.0, -40.0, -40.0]
    ax = Axis(Figure()[1, 1])
    AN.Plots.plotlayers!(ax, layers, depths)
    y = only(filter(p -> p isa HLines, ax.scene.plots))[1][]
    @test y ≈ [-90.0, -50.0] # half a row below each layer's first channel

    ax = Axis(Figure()[1, 1]) # a regular range keeps its own step
    AN.Plots.plotlayers!(ax, ["L1", "L1", "L2/3", "L2/3"], -100.0:20.0:-40.0)
    @test only(filter(p -> p isa HLines, ax.scene.plots))[1][] ≈ [-70.0] # between -80 and -60

    ax = Axis(Figure()[1, 1]) # descending depths, as a vector, keep the range's sign
    AN.Plots.plotlayers!(ax, ["L1", "L2/3", "L2/3"], [-40.0, -60.0, -80.0])
    @test only(filter(p -> p isa HLines, ax.scene.plots))[1][] ≈ [-50.0]
end

@testset "structure meshes" begin
    n = 21
    ball = [Float64(sum(abs2, (i, j, k) .- 11) ≤ 36) for i in 1:n, j in 1:n, k in 1:n] # radius 6 voxels
    m = AN.Plots.structuremesh(ball, [1:n, 1:n, .-(1:n)]) # axes as plotreferencevolume sets them
    @test !isempty(AN.Plots.GeometryBasics.faces(m))
    @test all(v -> 5 < hypot((v .- (11, 11, -11))...) < 7, AN.Plots.GeometryBasics.coordinates(m))
end
