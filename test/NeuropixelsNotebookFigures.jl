import AllenNeuropixels as AN
using WGLMakie
using DataFrames
using Statistics
using FileIO

session_table = AN.getsessiontable()

function numprobes(sessionid)
	session_table[session_table[!, :id].==sessionid, :probe_count][1]
end

function targetintersections(sessionid)
	s = session_table[session_table.id.==sessionid, :ecephys_structure_acronyms][1]
	targets = ["VISp", "VISl", "VISrl", "VISal", "VISpm", "VISam", "CA1", "CA3", "DG", "SUB", "ProS", "LGd", "LP", "APN"];
	structures = eachmatch(r"'(.*?)'", s)
	structures = [replace(x.match, r"'"=>s"") for x ∈ structures]
	targetsintersected = [x ∈ structures for x ∈ targets]
	return mean(targetsintersected)
end;

metrics = AN.getunitanalysismetricsbysessiontype("brain_observatory_1.1")

session_metrics = combine(groupby(metrics, :ecephys_session_id),
		:ecephys_session_id => numprobes∘unique => :num_probes,
		:ecephys_session_id => targetintersections∘unique => :target_intersections,
		:has_lfp_data=>all,
		:genotype => (x->all(isequal.(x, ("wt/wt",)))) => :is_wt,
		:max_drift=>median∘skipmissing,
		:d_prime=>median∘skipmissing,
		:isolation_distance=>median∘skipmissing,
		:silhouette_score=>median∘skipmissing,
		:snr=>median∘skipmissing)

oursession = subset(session_metrics,
						:num_probes => ByRow(==(6)),
						:target_intersections => ByRow(>(0.9)),
						:is_wt => ByRow(==(true)),
						:has_lfp_data_all => ByRow(==(true)),
						:max_drift_median_skipmissing => ByRow(<(25.0)))

session = AN.Session(oursession.ecephys_session_id[1])




# --------------------- Do the plots and save to dropbox --------------------- #
AN.Plots.formattedreferencevolume(session, "C:\\Users\\Brendan\\Dropbox (Sydney Uni Student)\\AllenNeuropixels\\probelocations.html")