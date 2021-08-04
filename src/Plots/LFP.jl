using DimensionalData
using Statistics

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
