function ecephyscache()
    ecephys_project_cache.EcephysProjectCache.from_warehouse(manifest=ecephysmanifest)
end

# These will take a while to download
function getsessiontable()
    @info "Please wait, this can take a few seconds"
    CSV.read(IOBuffer(ecephyscache().get_session_table().to_csv()), DataFrame);
end
export getsessiontable

function getprobes()
    CSV.read(IOBuffer(ecephyscache().get_probes().to_csv()), DataFrame);
end
export getprobes

function getchannels()
    CSV.read(IOBuffer(ecephyscache().get_channels().to_csv()), DataFrame);
end
export getchannels


function getunits(; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5)
    str = ecephyscache().getunits(filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum).to_csv()
    CSV.read(IOBuffer(str), DataFrame);
end
export getunits


function getunitanalysismetricsbysessiontype(session_type; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5) # Yeah thats python
    str = ecephyscache().get_unit_analysis_metrics_by_session_type(session_type,
                            filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum).to_csv()
    CSV.read(IOBuffer(str), DataFrame);
end
export getunitanalysismetricsbysessiontype


function getallunitmetrics() # This one is really slow
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
export getallunitmetrics


function getsessiondata(session_id::Int; filter_by_validity=true, amplitude_cutoff_maximum = 0.1, presence_ratio_minimum = 0.9, isi_violations_maximum = 0.5)
    ecephyscache().get_session_data(session_id; filter_by_validity=filter_by_validity,
                            amplitude_cutoff_maximum=amplitude_cutoff_maximum,
                            presence_ratio_minimum=presence_ratio_minimum,
                            isi_violations_maximum=isi_violations_maximum)
end
export getsessiondata


abstract type AbstractSession end

struct Session <: AbstractSession
    pyObject
end
export Session
Session(session_id::Int; kwargs...) = Session(getsessiondata(session_id; kwargs...))
getid(S::AbstractSession) = S.pyObject.ecephys_session_id
getprobes(S::AbstractSession) = CSV.read(IOBuffer(S.pyObject.probes.to_csv()), DataFrame)
getprobeids(S::AbstractSession) = getprobes(S)[!, :id]
getchannels(S::AbstractSession) = CSV.read(IOBuffer(S.pyObject.channels.to_csv()), DataFrame)
function getprobecoordinates(S::AbstractSession)
    c = subset(getchannels(S),              :anterior_posterior_ccf_coordinate => ByRow(!ismissing),
                                            :dorsal_ventral_ccf_coordinate => ByRow(!ismissing),
                                            :left_right_ccf_coordinate => ByRow(!ismissing))
    x = c[!, :anterior_posterior_ccf_coordinate]
    y = c[!, :dorsal_ventral_ccf_coordinate]
    z = c[!, :left_right_ccf_coordinate]
    return (x, y, z)
end

