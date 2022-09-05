
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
