import AllenNeuropixels as AN
using DataFrames

𝒮 = AN.get_session_table()
sessionid = 𝒮[1, :id] # Chelsea: 732592105

S = AN.Session(sessionid)
probes = AN.getprobes(S)
probeid = probes[1, :id]
data = AN.getlfp(S, probeid)
channels = subset(AN.get_channels(),    :ecephys_probe_id           => ByRow(==(probeid)),
                                        :ecephys_structure_id       => ByRow(!ismissing),
                                        :ecephys_structure_acronym  => ByRow(x->x!=="grey"))

AN.plotreferencevolume(S)