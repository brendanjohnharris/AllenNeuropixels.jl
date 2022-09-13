
function spikeraster(times::AbstractDict)

end


function spikeraster!(ax::Axis, spikes::AN.SpikeMatrix; markersize=15, marker=:vline, kwargs...)
    spy!(ax, dims(spikes, Ti)|>collect, 1:size(spikes, 2), spikes|>Array; markersize, marker)
    ax.xlabel = "Time(s)"
    ax.ylabel = "Unit"
    return ax
end

function spikeraster!(ax::Axis, spikes, stimuli::DataFrame; kwargs...)
    spikeraster!(ax, spikes; kwargs...)
    stimulusvlines!(ax, spikes, stimuli; stimcolor=nothing)
    return ax
end

function Makie.spy!(ax, x::AbstractVector{<:AbstractVector}; colormap=nothing, kwargs...)
    # Could color height by correlation
    if isnothing(colormap)
        [scatter!(ax, _x, zeros(length(_x)).+i; kwargs...) for (i, _x) in (enumerate)(x)]
    else
        [scatter!(ax, _x, zeros(length(_x)).+i; color=zeros(length(_x)).+i, colorrange=(1, length(x)), colormap, kwargs...) for (i, _x) in (enumerate)(x)]
    end
end
# function Makie.spy!(ax, x::Dict; kwargs...)
#     [scatter!(ax, _x, zeros(length(_x)).+i; kwargs...) for (i, _x) in (enumerate∘values)(x)]
# end


function plotspikebursts!(ax, X::AN.LFPVector, Y::Dict; marker=:vline, linecolor=:crimson, colormap=nothing, kwargs...)
    sims = zeros(length(Y))
    ts = Interval(extrema(dims(X, Ti))...)
    Y = deepcopy(Y)
    [(Y[i] = Y[i][Y[i] .∈ (ts,)]) for i in eachindex(Y)]
    # Sort by mean, standardised X value at each T
    X = (X .- mean(X))./std(X)
    for (i, k) in enumerate(eachindex(Y))
        sims[i] = mean(X[Ti(Near(Y[k]))])
    end
    idxs = sortperm(sims)
    times = collect(values(Y))[idxs]
    times = times[length.(times) .> 0]
    times = Vector{Vector}(times)
    spy!(ax, times; colormap, marker)
    lines!(ax, dims(X, Ti)|>collect, collect(length(times).*(X.-minimum(X))./(maximum(X) - minimum(X))); color=linecolor)
    ax.xlabel = "Time (s)"
    ax.ylabel = "Unit"
end
