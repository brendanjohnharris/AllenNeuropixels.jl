
function spikeraster(times::AbstractDict)

end


function spikeraster!(ax::Axis, spikes::AN.SpikeMatrix; kwargs...)
    spy(dims(spikes, Ti)|>collect, 1:size(spikes, 2), spikes|>Array, markersize=4)
    ax.xlabel = "Time(s)"
    ax.ylabel = "Unit"
    return ax
end

function spikeraster!(ax::Axis, spikes, stimuli::DataFrame; kwargs...)
    spikeraster!(ax, spikes; kwargs...)
    stimulusvlines!(ax, spikes, stimuli; stimcolor=nothing)
    return ax
end
