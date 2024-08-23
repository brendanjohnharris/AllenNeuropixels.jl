using DimensionalData
import DimensionalData as DD
using Statistics
using DSP
using StatsBase
using IntervalSets
import AllenNeuropixelsBase as ANB

function slidingcarpet(x, y, Z; size = (800, 800), kwargs...)
    # x is a vector of times labelling columns of Z, , y is a vector labelling rows of Z and Z is the matrix of time series data
    fig = Figure(; size = size)
    # sl_a = Makie.Slider(fig[1, 1], range = 0:0.01:10, startvalue = 3)
    # sl_t = Makie.Slider(fig[2, 1], range = 0:0.01:10, startvalue = 3)
    sla = labelslider!(fig, "Window", (2:length(x)))
    fig[1, 2] = sla.layout
    slb = labelslider!(fig, "Time     ", 0.01:0.005:1)
    fig[2, 2] = slb.layout

    u = lift(sla.slider.value, slb.slider.value) do sl...
        r = floor(length(x) * sl[1] ./ length(x)) # The number of time points to show
        s = ceil((length(x) - r) * sl[2]) # The starting position

        extrema([s, s + r - 1])
    end
    limits = lift(x -> (to_value(x)[1], to_value(x)[2], nothing, nothing), u)
    colsize!(fig.layout, 1, Relative(1 / 6))
    ax = Axis(fig[3, 2], limits = limits)
    heatmap!(ax, Z; kwargs...)
    xx = round.(LinRange(1, round(length(x), sigdigits = 1), 5), sigdigits = 1)
    xxx = round.(x[Int.(xx)], sigdigits = 2)
    xxx[ismissing.(xxx)] .= NaN
    #ax.xticks = (xx, string.(xxx))
    y[ismissing.(y)] .= ""
    ax.yticks = (1:size(Z, 2), y)
    return fig
end
export slidingcarpet

function slidingcarpet(X::DimArray; size = (800, 400), kwargs...)
    x = dims(X, 𝑡).val
    y = dims(X, Chan).val
    Z = Array(X)
    slidingcarpet(x, y, Z; size = size, kwargs...)
end

function neuroslidingcarpet(X::DimArray; size = (800, 400), kwargs...)
    time = dims(X, 𝑡).val
    channels = AN.getstructureacronyms(Meta.parse.(string.(dims(X, Chan).val)))
    fig = slidingcarpet(time, channels, X .- mean(Array(X), dims = 1);
                        size = size, kwargs...)
end

"""
We'll plot the columns of this array as traces of their mean values. `x` is the lag index of each row, `tracez` is the value with which to color each trace (column) and `X` is an array of trace data.
"""
function traces(x, X::AbstractArray; smooth = 10, clabel = nothing, tracez = nothing,
                colormap = cgrad([:cornflowerblue, :black], alpha = 0.2), domean = false,
                doaxis = true, linewidth = 1, alpha = 0.32, kwargs...) # ! Plot recipe?
    function MA(g)
        [i < smooth ? mean(g[begin:i]) : mean(g[(i - smooth + 1):i]) for i in 1:length(g)]
    end
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1]; kwargs...)
    y = StatsBase.median(X, dims = 2)[:]
    streamlines = vcat(mapslices(MA, X, dims = 1), fill(NaN, (1, size(X, 2))))[:]
    xx = repeat(x, outer = [1, size(X, 2)])
    xx = vcat(xx, fill(NaN, (1, size(xx, 2))))[:]
    if isnothing(tracez)
        colors = (:gray, alpha)
    else
        colors = repeat(tracez[:]', outer = [size(X, 1), 1])
        colors = Float64.(vcat(colors, fill(NaN, (1, size(colors, 2))))[:])
    end
    if doaxis
        vlines!(ax, [0], color = :gray42, linewidth = 2)
        hlines!(ax, [0], color = :gray42, linewidth = 2)
    end
    s = lines!(ax, xx, streamlines; color = colors, linewidth, colormap)
    if !isnothing(tracez)
        Colorbar(fig[1, 2], limits = extrema(colors[.!isnan.(colors)]),
                 colormap = cgrad(colormap, alpha = 0.7), flipaxis = true, label = clabel)
    end
    if domean
        lines!(ax, x, MA(y), linewidth = 5, color = :gray32)
    end
    return fig
end
export traces

function powerlawfit(_psd::AN.PSDMatrix)
    y = log10.(mean(_psd, dims = Chan))
    x = dims(_psd, 𝑓) |> collect
    x = repeat(log10.(x), 1, size(y, 2))
    y = Array(y)[:] # Flatten for regression
    x = x[:]
    X = [ones(size(y, 1)) x]
    b = X \ y
    c, r = b
    f = x -> (10^c) .* x .^ r
    return c, r, f
end

function plotLFPspectra(LFP::AbstractDimArray; slope = nothing,
                        position = Point2f([5, 1e-5]), fs = nothing, N = 500,
                        slopecolor = :crimson, kwargs...)
    times = collect(dims(LFP, 𝑡))
    if isnothing(fs)
        Δt = times[2] - times[1]
        all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
    else
        Δt = 1 / fs
    end
    fp = x -> welch_pgram(x, div(length(x), N), div(div(length(x), N), 2); fs = 1 / Δt,
                          window = nothing)
    P = [fp(Array(x)) for x in eachcol(LFP)]
    frq = P[1].freq # Should be pretty much the same for all columns?
    psd = hcat([p.power for p in P]...)
    psd = psd ./ (sum(psd, dims = 1) .* (frq[2] - frq[1]))
    psd = DimArray(psd, (𝑓(frq), dims(LFP, Chan)))
    fig = traces(frq, Array(psd); xlabel = "frq (Hz)", ylabel = "Ŝ",
                 title = "Normalised power spectral density", smooth = 1,
                 yscale = Makie.log10, xscale = Makie.log10, doaxis = false, domean = false,
                 yminorgridvisible = false, kwargs...)
    if !isnothing(slope)
        _psd = psd[𝑓(DD.Between(slope...))]
        c, r, f = powerlawfit(_psd)
        lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color = slopecolor,
               linewidth = 5)
        text!(L"$\alpha$= %$(round(r, sigdigits=2))", position = Point2f0(position),
              fontsize = 40)
    end
    return fig
end

function plotLFPspectra(session, probeid, LFP::AbstractDimArray; slope = nothing,
                        position = Point2f([5, 1e-5]), slopecolor = :crimson, kwargs...)
    # Calculate the power spectrum of each column of the LFP array
    times = collect(dims(LFP, 𝑡))
    Δt = times[2] - times[1]
    all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
    # frq = rfftfreq(length(times), 1/Δt)
    # 𝑍 = rfft(Array(LFP), 1)
    # A = abs.(𝑍)
    # psd = (Δt/length(times))*A.^2
    fp = x -> welch_pgram(x, div(length(x), 500), div(div(length(x), 500), 2); fs = 1 / Δt,
                          window = nothing)
    P = [fp(Array(x)) for x in eachcol(LFP)]
    frq = P[1].freq # Should be pretty much the same for all columns?
    psd = hcat([p.power for p in P]...)
    psd = psd ./ (sum(psd, dims = 1) .* (frq[2] - frq[1]))
    psd = DimArray(psd, (𝑓(frq), dims(LFP, Chan)))
    depths = ANB.getchanneldepths(session, probeid, collect(dims(psd, Chan)))
    fig = traces(frq, Array(psd); tracez = -depths, xlabel = "frq (Hz)", ylabel = "Ŝ",
                 title = "Normalised power spectral density", clabel = "Depth", smooth = 1,
                 yscale = Makie.log10, doaxis = false, domean = false,
                 yminorgridvisible = false, kwargs...)
    if !isnothing(slope)
        _psd = psd[𝑓(DD.Between(slope...))]
        c, r, f = powerlawfit(_psd)
        lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color = slopecolor,
               linewidth = 5)
        text!(L"$\alpha$= %$(round(r, sigdigits=2))", position = Point2f0(position),
              fontsize = 40)
    end
    return fig
end

function plotLFPspectra(session, probeids::Vector, LFP::Vector; slopecolor = nothing,
                        kwargs...)
    @assert length(probeids) == length(LFP)
    times = collect(dims(LFP[1], 𝑡))
    Δt = times[2] - times[1]
    @assert all(Δt .≈ diff(times))
    fp = x -> welch_pgram(x, div(length(x), 200), div(div(length(x), 200), 2); fs = 1 / Δt,
                          window = nothing)
    A = Vector(undef, length(LFP))
    P = fp(Array(LFP[1][:, 1]))
    frq = P[1].freq # Should be pretty much the same for all channels?
    @time for x in 1:length(LFP)
        P = [fp(Array(x)) for x in eachcol(LFP)]
        psd = hcat([p.power for p in P]...)
        psd = psd ./ (sum(psd, dims = 1) .* (frq[2] - frq[1]))
        A[x] = DimArray(psd, (𝑓(frq), dims(LFP, Chan)))
    end
    depths = [ANB.getchanneldepths(session, probeids[x], collect(dims(A[x], Chan)))
              for x in 1:length(probeids)]
    A = hcat(A...)
    depths = vcat(depths...)
    # idxs = .!isnan.(depths)
    # AC = AC[:, idxs]
    # depths = depths[idxs]
    traces(timelags, AC; tracez = -depths, xlabel = "𝜏 (s)", ylabel = "𝜌(𝜏)",
           title = "Autocorrelation of LFP timeseries", clabel = "Depth", smooth = 10,
           domean = false, colormap = cgrad([:cornflowerblue, :black], alpha = 0.1),
           kwargs...)
end
export plotLFPspectra

function stackedtraces!(ax::Axis, X::AN.LFPMatrix; offset = 0.5, stimulus = nothing,
                        stimcolor = nothing, rev = false, kwargs...)
    inc = 0
    rev && (X = -X)
    #c = 0.0#0.1*mean(diff.(extrema.(eachcol(X)).|>collect))[1]
    # c = offset*mean(diff.(extrema.(eachslice(X, dim=Chan)).|>collect))[1]
    #data = [(inc += (c + minimum(X[:, i] .- X[:, i+1])[1]); X[:, i] .+ inc) for i ∈ 1:size(X, 2)-1]
    X = X .- minimum.(eachcol(X))'
    channels = dims(X, Chan)

    if isnothing(offset)
        X = X ./ (maximum.(eachcol(X))' .- minimum.(eachcol(X))')
        c = 1:size(X, Chan)
    else
        c = zeros(size(X, Chan))
        for i in 2:size(X, Chan)
            c[i] = c[i - 1] + maximum(X[:, i - 1] .- X[:, i]) + offset * mean(X)
        end
    end
    data = [X[:, i] .+ c[i] for i in 1:size(X, 2)]
    Makie.lines!.((ax,), (dims(X, 𝑡) |> collect,), data |> collect; kwargs...)
    # hlines!(ax, c)
    ax.yticks = (mean.(data), string.(channels))
    ax.yticklabelrotation = 0 # π/2
    ax.xlabel = "time (s)"
    ax.ylabel = "channel"
    ax.yticklabelsvisible = true
    ax.yticksvisible = false
    if stimulus isa DataFrame
        stimulusvlines!(ax, X, stimulus; stimcolor)
    end
    rev && (ax.yreversed = true)
end

function stackedtraces(X::AN.LFPMatrix; kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1])
    stackedtraces!(ax, X; kwargs...)
    display(fig)
    (fig, ax)
end

function stimulusvlines!(ax, X, stimulus; stimcolor = nothing)
    times = Interval(extrema(dims(X, 𝑡))...)
    starts = stimulus.start_time
    stops = stimulus.stop_time
    startidxs = starts .∈ (times,)
    starts = starts[startidxs]
    stops = stops[stops .∈ (times,)]
    if stimcolor isa Symbol
        linez = stimulus[startidxs, stimcolor]
        if eltype(linez) <: AbstractString
            linez = Vector{Any}(Meta.parse.(linez))
            linez[Symbol.(string.(linez)) .== :null] .= 0.0
            linez = Vector{Float32}(linez)
        end
        Makie.vlines!(ax, starts, color = linez, linewidth = 5)
    else
        Makie.vlines!(ax, starts, color = :green, linewidth = 5)
        Makie.vlines!(ax, stops, color = :red, linewidth = 5)
    end
end

function regiontraces(X::AN.LFPMatrix; kwargs...)
    channels = dims(X, Chan) |> collect
    structures = AN.getstructureacronyms(channels)
    fig, ax = stackedtraces(X; kwargs...)
    ax.yticks = (ax.yticks.val[1], structures)
end

function plotwaveletvariability(res::AN.LogWaveletMatrix, L::Function; downsample = 10)
    x = dims(res, Log𝑓) |> collect
    fig = traces((x), collect(res[1:downsample:end, :])'; domean = true, alpha = 0.01) # Plots the wavelet spectrum at each time against the mean wavelet spectrum
    ax = current_axis()
    lines!(ax, (x), L.(x), color = :crimson, linewidth = 2)
    ax.xlabel = "log10(frequency)"
    return fig
end

function Makie.heatmap!(ax, B::AN.Burst; kwargs...)
    x = dims(B.mask, 𝑡) |> collect
    y = dims(B.mask, Log𝑓) |> collect
    heatmap!(ax, x, y, B.mask |> Array; kwargs...)
end
function Makie.heatmap(B::AN.Burst; kwargs...)
    (ax = Axis(Figure()[1, 1]); heatmap!(ax, B; kwargs...); current_figure())
end

function plotfit!(ax, b::AN.Burst; kwargs...)
    x = dims(b.mask, 𝑡) |> collect
    y = dims(b.mask, Log𝑓) |> collect
    surface!(ax, x, y .|> exp10, b.mask |> Array; kwargs...)
    model = xy -> AN.gaussian2(xy, b.fit.param)
    x = LinRange(extrema(x)..., 10)
    y = LinRange(extrema(y)..., 10)
    ax.xlabel = "Time (s)"
    ax.ylabel = "Frequency (Hz)"
    wireframe!(ax, x, y .|> exp10, model.(Iterators.product(x, y)))
end
function plotfit(B::AN.Burst; kwargs...)
    (ax = Axis3(Figure()[1, 1]); plotfit!(ax, B; kwargs...); current_figure())
end

function plotfit!(ax, res::AN.LogWaveletMatrix, B::AN.BurstVector;
                  N = min(10000, size(res, 1)), downsample = 1, colorbar = 1,
                  colorrange = extrema(res[1:downsample:N, :]), contourcolor = :crimson,
                  strokecolor = :gray32,
                  kwargs...)
    ctitle = "ΔS (a.u.)"
    ax.xlabel = "Time (s)"
    ax.ylabel = "Frequency (Hz)"#, yscale=Makie.pseudolog10);
    t, freqs, res = decompose(res)
    res = Array(res)
    ts = ClosedInterval(extrema(t[1:downsample:N])...)
    fs = ClosedInterval(extrema(exp10.(freqs))...)
    p = Makie.heatmap!(ax, t[1:downsample:N], exp10.(freqs), res[1:downsample:N, :];
                       colorrange, kwargs...)
    colorbar > 0 && Colorbar(current_figure()[colorbar, 2], p; label = ctitle)

    for b in B[AN.peaktime.(B) .∈ (ts,)]
        t, freqs, res = decompose(b.mask)
        res = Array(res)
        contour!(ax, t[1:downsample:end], exp10.(freqs), res[1:downsample:end, :];
                 color = contourcolor, linewidth = 2)
        vs = collect(Iterators.product((b |> AN.mask |> AN.waveletmatrix |> dims .|>
                                        extrema)...))
        poly!(Point2f.(vs[[1, 2, 4, 3]]); color = RGBA(0, 0, 0, 0), strokecolor,
              strokewidth = 4)
    end
    ax.limits = (extrema(ts), extrema(fs))
    return ax
end
function plotfit(res::AN.LogWaveletMatrix, B::AN.BurstVector; kwargs...)
    (ax = Axis(Figure()[1, 1]); plotfit!(ax, res, B; kwargs...); current_figure())
end

function _plotfit!(ax::Axis, b::AN.Burst; downsample = 1,
                   contourcolor = (:cornflowerblue, 0.4), boxcolor = :crimson,
                   strokewidth = 4)
    t, freqs, res = decompose(b.mask)
    contour!(ax, t[1:downsample:end], freqs, res[1:downsample:end, :]; color = contourcolor,
             linewidth = 2)
    vs = collect(Iterators.product((b |> AN.mask |> dims .|> extrema)...))
    poly!(ax, Point2f.(vs[[1, 2, 4, 3]]); color = Makie.RGBA(0, 0, 0, 0),
          strokecolor = boxcolor, strokewidth)
end

function compareburststatistics(B::AN.BurstVector, Bs::AN.BurstVector, f::Function,
                                label1 = "Real", label2 = "Surrogate"; boundary = nothing,
                                kwargs...)
    fig = Figure()
    ax = Axis(fig[1, 1], ylabel = "Probability density", xlabel = string(f))
    ax2 = Axis(fig[2, 1], ylabel = "Count", xlabel = "log10(p-value)")
    d = f.(B)
    ds = f.(Bs)
    kwargs = ()
    if !isnothing(boundary)
        kwargs = (; boundary = boundary)
    end

    p1 = hist!(ax, d; label = label1, kwargs...)
    p2 = hist!(ax, ds; label = label2, kwargs...)
    axislegend(ax)

    ps = AN.pvalues(B, Bs, f)
    p3 = hist!(ax2, log10.(ps .+ eps()))
    display(fig)
    return fig, ax, (p1, p2, p3)
end
