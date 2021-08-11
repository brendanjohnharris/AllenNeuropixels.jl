using IntervalSets

function brainobservatorycache(manifest_file=brainobservatorymanifest)
   brain_observatory_cache.BrainObservatoryCache(manifest_file=manifest_file)
end

function getophysexperiments(; kwargs...)
    brainobservatorycache().get_ophys_experiments(; kwargs...)
end


function getspatialgrating(; height=100, aspect_ratio=2, ori=45, pix_per_cycle=1/0.04, phase=0, p2p_amp=0.8, baseline=p2p_amp/2)
    if any([ori, pix_per_cycle, phase, p2p_amp, baseline] .== ["null"]) # They showed blank (gray) gratings
        f = stimulus_info.get_spatial_grating(; height, aspect_ratio, ori=45, pix_per_cycle=height, phase=0, p2p_amp=0.8, baseline=0)
        f .= 0.4 # A gray color, mean of sine between 0 and 0.8
    end
    stimulus_info.get_spatial_grating(; height, aspect_ratio, ori, pix_per_cycle, phase, p2p_amp, baseline)
end

function getdriftinggrating(t::Real; temporal_frequency=2, kwargs...)
    d = kwargs |> Dict{Symbol,Any}
    d[:phase] =  pop!(d, :phase, 0) + (t*temporal_frequency)%1 # ! This is how they do it in the allensdk
    getspatialgrating(; d...)
end

function getnaturalmovie(stimulus)
    boc = getophysexperiments(stimuli=[stimulus])[1] # Just need one of these sessions to get the template
    ds = brainobservatorycache().get_ophys_experiment_data(boc["id"])
    template = ds.get_stimulus_template(stimulus)
end

function getnaturalstimulusframes(session, times)
    df = getstimuli(session, times)
    movies = unique(df.stimulus_name)
    @assert all(.∈(movies , (["natural_movie_one", "natural_movie_three", "natural_scenes"],)))
    movieframes = [getnaturalmovie(movie) for movie in movies]
    frames = Vector{Matrix{Float32}}(undef, length(times))
    for t in 1:length(times)
        movieidx = indexin([df.stimulus_name[t]], movies)[1]
        frame = Int(Meta.parse(df.frame[t])+1) # Goddamn python frames start at 0
        if frame < 1
            frames[t] = fill(1.0, size(frames[t-1]))
        else
            frames[t] = movieframes[movieidx][frame, :, :]./256
        end
    end
    return DimArray(frames, (Dim{:time}(times)))
end

getnaturalmovieframes = getnaturalstimulusframes

function getdriftinggratingframes(session, times)
    df = getstimuli(session, times)
    height = 100
    aspect_ratio = 2
    frames = [Matrix{Float32}(undef, height, aspect_ratio*height) for t ∈ 1:length(times)]
    for t ∈ 1:lastindex(times)
        if any((df.temporal_frequency[t], df.orientation[t], df.contrast) .== ["null"])
            frames[t] .= 0.4
        else
            frames[t] = getdriftinggrating(times[t] - df.start_time[t];
                    temporal_frequency=df.temporal_frequency[t]|>Meta.parse,
                    height,
                    aspect_ratio,
                    ori=df.orientation[t]|>Meta.parse,
                    pix_per_cycle=1/Meta.parse(df.spatial_frequency[t]),
                    p2p_amp=df.contrast[t]|>Meta.parse)
        end
    end
    return DimArray(frames, (Dim{:time}(times)))
end

function getstaticgratingframes(session, times, spontaneous=false)
    df = getstimuli(session, times)
    height = 100
    aspect_ratio = 2
    frames = [Matrix{Float32}(undef, height, aspect_ratio*height) for t ∈ 1:length(times)]
    for t ∈ 1:lastindex(times)
        if any((df.orientation[t], df.contrast) .== ["null"]) || spontaneous
            frames[t] .= 0.4
        else
            frames[t] = getspatialgrating(;
                    height,
                    aspect_ratio,
                    ori=df.orientation[t]|>Meta.parse,
                    pix_per_cycle=1/Meta.parse(df.spatial_frequency[t]),
                    p2p_amp=df.contrast[t]|>Meta.parse)
        end
    end
    return DimArray(frames, (Dim{:time}(times)))
end

function _getstimulusframes(session, times, stimulus)
    if stimulus ∈ ["static_gratings", "spontaneous"] # Spontaneous is just a mean-luminance blank screen]
        getstaticgratingframes(session, times, stimulus=="spontaneous")
    elseif stimulus == "drifting_gratings"
        getdriftinggratingframes(session, times)
    elseif stimulus ∈ ["natural_movie_one", "natural_movie_three", "natural_scenes"]
        getnaturalstimulusframes(session, times)
    end
end

function getstimulusframes(session, times)
    epochs = getepochs(session)
    epochs = epochs[[any(times .∈ x.start_time..x.stop_time) for x ∈ eachrow(epochs)], :] # Get epochs that contain the times
    epochtimes = [[t for t ∈ times if t ∈ x.start_time..x.stop_time] for x ∈ eachrow(epochs)]
    frames = [_getstimulusframes(session, epochtimes[x], epochs.stimulus_name[x]) for x ∈ 1:length(epochtimes)]
    frames = DimArray(vcat(frames...), (Dim{:time}(times)))
end
