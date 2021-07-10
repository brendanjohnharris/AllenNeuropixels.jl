function brainobservatorycache(manifest_file=brainobservatorymanifest)
   brain_observatory_cache.BrainObservatoryCache(manifest_file=manifest_file)
end

function getophysexperiments(; kwargs...)
    brainobservatorycache().get_ophys_experiments(; kwargs...)
end


function getspatialgrating(; height=100, aspect_ratio=2, ori=45, pix_per_cycle=10, phase=0, p2p_amp=2, baseline=1)
    stimulus_info.get_spatial_grating(; height, aspect_ratio, ori, pix_per_cycle, phase, p2p_amp, baseline)
end

