function cache()
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=manifestdir)
end

# These will take a while to download
function get_session_table()
    types = [String, Int64, String, Float64, Symbol, String, Int, Int, Int, Vector]
    df = convertdataframe(cache().get_session_table(), types)
    return df
end
export get_session_table

function get_probes()
    types = [Int64, Float64, String, String, Float64, Bool, Int64, Int64, Vector]
    convertdataframe(cache().get_probes(), types)
end
export get_probes

function get_channels()
    types = [Int, Int, Int, Int, Int, Int, Int, Int, Symbol, Int, Float64, Symbol, Float64, Bool, Int64]
    convertdataframe(cache().get_channels(), types)
end
export get_channels

# function get_units()
#     types = [Int, Int, Int, Int, Int, Int, Int, Int, Symbol, Int, Float64, Symbol, Float64, Bool, Int64]
#     convertdataframe(cache().get_units(), types)
# end
# export get_units
