function cache()
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=manifestpath)
end

# These will take a while to download
function get_session_table(file="sessions.csv", dir=datadir) # Really, this should be called with no arguments
    if !isfile(abspath(dir, file))
        cache().get_sessions() # Download data
    end
    loaddataframe(file, dir)
end
export get_session_table

function get_probes(file="probes.csv", dir=datadir)
    if !isfile(abspath(dir, file))
        cache().get_probes()
    end
    loaddataframe(file, dir)
end
export get_probes

function get_channels(file="channels.csv", dir=datadir)
    if !isfile(abspath(dir, file))
        cache().get_channels()
    end
    loaddataframe(file, dir)
end
export get_channels


function get_units(file="units.csv", dir=datadir)
    if !isfile(abspath(dir, file))
        cache().get_units()
    end
    loaddataframe(file, dir)
end
export get_units

