import AllenNeuropixels as AN
using DataFrames

sessions = AN.get_session_table()

session = subset(sessions, :sex => ByRow(==("M")), :full_genotype => ByRow(x->contains(x, "Sst")), :session_type => ByRow(==("brain_observatory_1.1")), :ecephys_structure_acronyms => ByRow(x->contains(x, "VISl")))

S = AN.Session(session[1, :id])
probes = AN.getprobeids(S)
data = AN.getlfp(S, probes[1])