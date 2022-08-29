using DimensionalData
import DimensionalData as DD
using Statistics
using DSP
using StatsBase
using IntervalSets

function slidingcarpet(x, y, Z; resolution=(800, 800), kwargs...)
    # x is a vector of times labelling columns of Z, , y is a vector labelling rows of Z and Z is the matrix of time series data
    fig = Figure(;resolution=resolution)
    # sl_a = Makie.Slider(fig[1, 1], range = 0:0.01:10, startvalue = 3)
    # sl_t = Makie.Slider(fig[2, 1], range = 0:0.01:10, startvalue = 3)
    sla = labelslider!(fig, "Window", (2:length(x)))
    fig[1, 2] = sla.layout
    slb = labelslider!(fig, "Time     ", 0.01:0.005:1)
    fig[2, 2] = slb.layout

    u = lift(sla.slider.value, slb.slider.value) do sl...
        r = floor(length(x)*sl[1]./length(x)) # The number of time points to show
        s = ceil((length(x) - r)*sl[2]) # The starting position

        extrema([s, s+r-1])
    end
    limits = lift(x -> (to_value(x)[1], to_value(x)[2], nothing, nothing), u)
    colsize!(fig.layout, 1, Relative(1/6))
    ax = Axis(fig[3, 2], limits=limits)
    heatmap!(ax, Z; kwargs...)
    xx = round.(LinRange(1, round(length(x), sigdigits=1), 5), sigdigits=1)
    xxx = round.(x[Int.(xx)], sigdigits=2)
    xxx[ismissing.(xxx)] .= NaN
    #ax.xticks = (xx, string.(xxx))
    y[ismissing.(y)] .= ""
    ax.yticks = (1:size(Z, 2), y)
    return fig
end
export slidingcarpet


function slidingcarpet(X::DimArray; resolution=(800, 400), kwargs...)
    x = dims(X, Ti).val
    y = dims(X, :channel).val
    Z = Array(X)
    slidingcarpet(x, y, Z; resolution=resolution, kwargs...)
end

function neuroslidingcarpet(X::DimArray; resolution=(800, 400), kwargs...)
    time = dims(X, Ti).val
    channels = AN.getstructureacronyms(Meta.parse.(string.(dims(X, :channel).val)))
    fig = slidingcarpet(time, channels, X.-mean(Array(X), dims=1); resolution=resolution, kwargs...)
end



"""
We'll plot the columns of this array as traces of their mean values. `x` is the lag index of each row, `tracez` is the value with which to color each trace (column) and `X` is an array of trace data.
"""
function traces(x, X::AbstractArray; smooth=10, clabel=nothing, tracez=nothing, colormap=cgrad([:cornflowerblue, :black], alpha=0.2), domean=false, doaxis=true, linewidth=1, kwargs...) # ! Plot recipe?
    MA(g) = [i < smooth ? mean(g[begin:i]) : mean(g[i-smooth+1:i]) for i in 1:length(g)]
    fig = Figure(resolution=(800, 400))
    ax = Axis(fig[1, 1]; kwargs...)
    y = StatsBase.median(X, dims=2)[:]
    streamlines = vcat(mapslices(MA, X, dims=1), fill(NaN, (1, size(X, 2))))[:]
    xx = repeat(x, outer=[1, size(X, 2)])
    xx = vcat(xx, fill(NaN, (1, size(xx, 2))))[:]
    if isnothing(tracez)
        colors = (:gray, 0.32)
    else
        colors = repeat(tracez[:]', outer=[size(X, 1), 1])
        colors = Float64.(vcat(colors, fill(NaN, (1, size(colors, 2))))[:])
    end
    if doaxis
        vlines!(ax, [0], color = :gray42, linewidth=2)
        hlines!(ax, [0], color = :gray42, linewidth=2)
    end
    s = lines!(ax, xx, streamlines; color=colors, linewidth, colormap)
    if !isnothing(tracez)
        Colorbar(fig[1, 2], limits=extrema(colors[.!isnan.(colors)]), colormap=cgrad(colormap, alpha=0.7), flipaxis = true, label=clabel)
    end
    if domean
        lines!(ax, x, MA(y), linewidth=5, color=:gray32)
    end
    return fig
end
export traces


function powerlawfit(_psd::AN.PSDMatrix)
    y = log10.(mean(_psd, dims=Dim{:channel}))
    x = dims(_psd, Dim{:𝑓}) |> collect
    x = repeat(log10.(x), 1, size(y, 2))
    y = Array(y)[:] # Flatten for regression
    x = x[:]
    X = [ones(size(y, 1)) x]
    b = X\y
    c, r = b
    f = x->(10^c).*x.^r
    return c, r, f
end


function plotLFPspectra(LFP::AbstractDimArray; slope=nothing, position=Point2f([5, 1e-5]), kwargs...)
    times = collect(dims(LFP, Ti))
    Δt = times[2] - times[1]
    all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
    fp = x -> welch_pgram(x, div(length(x), 1000), div(div(length(x), 1000), 2); fs=1/Δt, window=nothing)
    P = [fp(Array(x)) for x ∈ eachcol(LFP)]
    𝑓 = P[1].freq # Should be pretty much the same for all columns?
    psd = hcat([p.power for p ∈ P]...)
    psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
    psd = DimArray(psd, (Dim{:𝑓}(𝑓), dims(LFP, :channel)))
    fig = traces(𝑓, Array(psd); xlabel="𝑓 (Hz)", ylabel="Ŝ", title="Normalised power spectral density", smooth=1, yscale=Makie.log10, doaxis=false, domean=false, yminorgridvisible=false, kwargs...)
    if !isnothing(slope)
        _psd = psd[Dim{:𝑓}(DD.Between(slope...))]
        c, r, f = powerlawfit(_psd)
        lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color=:red, linewidth=5)
        text!(L"$\alpha$= %$(round(r, sigdigits=2))", position=Point2f0(position), textsize=40)
    end
    return fig
end


function plotLFPspectra(session, probeid, LFP::AbstractDimArray; slope=nothing, position=Point2f([5, 1e-5]), kwargs...)
    # Calculate the power spectrum of each column of the LFP array
    times = collect(dims(LFP, Ti))
    Δt = times[2] - times[1]
    all(Δt .≈ diff(times)) || @warn "Violated assumption: all(Δt .≈ diff(times))"
    # 𝑓 = rfftfreq(length(times), 1/Δt)
    # 𝑍 = rfft(Array(LFP), 1)
    # A = abs.(𝑍)
    # psd = (Δt/length(times))*A.^2
    fp = x -> welch_pgram(x, div(length(x), 1000), div(div(length(x), 1000), 2); fs=1/Δt, window=nothing)
    P = [fp(Array(x)) for x ∈ eachcol(LFP)]
    𝑓 = P[1].freq # Should be pretty much the same for all columns?
    psd = hcat([p.power for p ∈ P]...)
    psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
    psd = DimArray(psd, (Dim{:𝑓}(𝑓), dims(LFP, :channel)))
    depths = AN.getchanneldepths(session, probeid, collect(dims(psd, :channel)))
    fig = traces(𝑓, Array(psd); tracez=depths, xlabel="𝑓 (Hz)", ylabel="Ŝ", title="Normalised power spectral density", clabel="Depth", smooth=1, yscale=Makie.log10, doaxis=false, domean=false, yminorgridvisible=false, kwargs...)
    if !isnothing(slope)
        _psd = psd[Dim{:𝑓}(DD.Between(slope...))]
        c, r, f = powerlawfit(_psd)
        lines!(LinRange(slope..., 100), f(LinRange(slope..., 100)), color=:red, linewidth=5)
        text!(L"$\alpha$= %$(round(r, sigdigits=2))", position=Point2f0(position), textsize=40)
    end
    return fig
end

function plotLFPspectra(session, probeids::Vector, LFP::Vector; kwargs...)
    @assert length(probeids) == length(LFP)
    times = collect(dims(LFP[1], Ti))
    Δt = times[2] - times[1]
    @assert all(Δt .≈ diff(times))
    fp = x -> welch_pgram(x, div(length(x), 1000), div(div(length(x), 1000), 2); fs=1/Δt, window=nothing)
    A = Vector(undef, length(LFP))
    P = fp(Array(LFP[1][:, 1]))
    𝑓 = P[1].freq # Should be pretty much the same for all channels?
    @time for x in 1:length(LFP)
        P = [fp(Array(x)) for x ∈ eachcol(LFP)]
        psd = hcat([p.power for p ∈ P]...)
        psd = psd./(sum(psd, dims=1).*(𝑓[2] - 𝑓[1]))
        A[x] = DimArray(psd, (Dim{:𝑓}(𝑓), dims(LFP, :channel)))
    end
    depths =[AN.getchanneldepths(session, probeids[x], collect(dims(A[x], :channel))) for x ∈ 1:length(probeids)]
    A = hcat(A...)
    depths = vcat(depths...)
    # idxs = .!isnan.(depths)
    # AC = AC[:, idxs]
    # depths = depths[idxs]
    traces(timelags, AC; tracez=depths, xlabel="𝜏 (s)", ylabel="𝜌(𝜏)", title="Autocorrelation of LFP timeseries", clabel="Depth", smooth=10, domean=false, colormap=cgrad([:cornflowerblue, :black], alpha=0.1), kwargs...)
end
export plotLFPspectra



function stackedtraces(X::AN.LFPMatrix; offset=0.5, stimulus=nothing, stimcolor=nothing, kwargs...)
    inc = 0
    #c = 0.0#0.1*mean(diff.(extrema.(eachcol(X)).|>collect))[1]
    # c = offset*mean(diff.(extrema.(eachslice(X, dim=Dim{:channel})).|>collect))[1]
    #data = [(inc += (c + minimum(X[:, i] .- X[:, i+1])[1]); X[:, i] .+ inc) for i ∈ 1:size(X, 2)-1]
    c = zeros(size(X, Dim{:channel}))
    for i in 2:size(X, Dim{:channel})
        y = X[:, i]
        x = X[:, i-1]
        c[i] = c[i-1] - maximum(x) + minimum(y)
    end
    c = c .*offset
    data = [X[:, i] .+ c[i] for i ∈ 1:size(X, 2)]
    fig = Figure(); ax = Axis(fig[1, 1], xlabel="time (s)", ylabel="channel", yticklabelsvisible=false, yticksvisible=false); Makie.lines!.((ax,), (dims(X, Ti)|>collect,), data|>collect; kwargs...)


    if stimulus isa DataFrame
        stimulusvlines!(ax, X, stimulus; stimcolor)
    end

    display(fig); (fig, ax)
end

function stimulusvlines!(ax, X, stimulus; stimcolor=nothing)
    times = Interval(extrema(dims(X, Ti))...)
    starts = stimulus.start_time
    stops = stimulus.stop_time
    startidxs = starts .∈ (times,)
    starts = starts[startidxs]
    stops = stops[stops .∈ (times,)]
    if stimcolor isa Symbol
        linez = stimulus[startidxs, stimcolor]
        if eltype(linez) <: AbstractString
            linez = Meta.parse.(linez)
            linez[linez .== :null] .= 0.0
            linez = Vector{Float32}(linez)
        end
        Makie.vlines!(ax, starts, color=linez, linewidth=5)
    else
        Makie.vlines!(ax, starts, color=:green, linewidth=5)
        Makie.vlines!(ax, stops, color=:red, linewidth=5)
    end
end
