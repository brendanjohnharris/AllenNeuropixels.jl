function cache()
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=manifestpath)
end

# These will take a while to download
function get_session_table()
    @info "Please wait, this can take a few seconds"
    CSV.read(IOBuffer(cache().get_session_table().to_csv()), DataFrame);
end
export get_session_table

function get_probes()
    CSV.read(IOBuffer(cache().get_probes().to_csv()), DataFrame);
end
export get_probes

function get_channels()
    CSV.read(IOBuffer(cache().get_channels().to_csv()), DataFrame);
end
export get_channels


function get_units(; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5)
    str = cache().get_units(filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum).to_csv()
    CSV.read(IOBuffer(str), DataFrame);
end
export get_units


function get_unit_analysis_metrics_by_session_type(session_type; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5) # Yeah thats python
    str = cache().get_unit_analysis_metrics_by_session_type(session_type,
                            filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum).to_csv()
    CSV.read(IOBuffer(str), DataFrame);
end
export get_unit_analysis_metrics_by_session_type


function get_all_unit_metrics() # This one is really slow
    metrics1 = get_unit_analysis_metrics_by_session_type("brain_observatory_1.1",
                            amplitude_cutoff_maximum = Inf,
                            presence_ratio_minimum = -Inf,
                            isi_violations_maximum = Inf)

    metrics2 = get_unit_analysis_metrics_by_session_type("functional_connectivity",
                            amplitude_cutoff_maximum = Inf,
                            presence_ratio_minimum = -Inf,
                            isi_violations_maximum = Inf)

    vcat(analysis_metrics1, analysis_metrics2)
end
export get_all_unit_metrics


function get_session_data(session_id::Int; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5)
    cache().get_session_data(session_id; filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum)
end
export get_session_data


abstract type AbstractSession end

struct Session <: AbstractSession
    pySession
end
export Session
Session(session_id::Int; kwargs...) = Session(get_session_data(session_id; kwargs...))
getid(S::AbstractSession) = S.pySession.ecephys_session_id
getprobes(S::AbstractSession) = CSV.read(IOBuffer(S.pySession.probes.to_csv()), DataFrame)
getprobeids(S::AbstractSession) = getprobes(S)[!, :id]
getlfp(S::AbstractSession, probeid::Int) = S.pySession.get_lfp(probeid)
