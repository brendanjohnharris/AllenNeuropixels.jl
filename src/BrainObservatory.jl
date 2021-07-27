function brainobservatorycache(manifest_file=brainobservatorymanifest)
   brain_observatory_cache.BrainObservatoryCache(manifest_file=manifest_file)
end

function getophysexperiments(; kwargs...)
    brainobservatorycache().get_ophys_experiments(; kwargs...)
end


function getspatialgrating(; height=100, aspect_ratio=2, ori=45, pix_per_cycle=10, phase=0, p2p_amp=2, baseline=1)
    stimulus_info.get_spatial_grating(; height, aspect_ratio, ori, pix_per_cycle, phase, p2p_amp, baseline)
end


function getnaturalmovie(stimulus="natural_movie_one")
    boc = getophysexperiments(stimuli=[stimulus])[1] # Just need one of these sessions to get the template
    ds = brainobservatorycache().get_ophys_experiment_data(boc["id"])
    template = ds.get_stimulus_template(stimulus)
end

function getnaturalmovieframes(session, times)
    df = getstimuli(session, times)
    movies = unique(df.stimulus_name)
    @assert all(.∈(movies , (["natural_movie_one", "natural_movie_three", "natural_scenes"],)))
    movieframes = [getnaturalmovie(movie) for movie in movies]
    res = size(movieframes[1])[2:3]
    frames = Array{typeof(movieframes[1][1])}(undef, length(times), res...)
    for t in 1:length(times)
        movieidx = indexin([df.stimulus_name[t]], movies)[1]
        frame = Int(Meta.parse(df.frame[t])+1) # Goddamn python frames start at 0
        if frame < 1
            frames[t, :, :] = fill(NaN, (1, size(frames)[2:3]...))
        else
            frames[t, :, :] = movieframes[movieidx][frame, :, :]
        end
    end
    return DimArray(frames, (Dim{:window}(times), DimensionalData.X(1:size(frames, 2)), DimensionalData.Y(1:size(frames, 3))))
end
