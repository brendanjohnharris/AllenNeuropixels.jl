import AllenNeuropixels as AN
using Profile

params = (;
          sessionid = 1044385384,
          stimulus = "spontaneous",
          structure = "VISl")

function proftest()
    session = AN.Session(params[:sessionid])
    probeids = AN.getprobestructures(session)
    LFP = AN.formatlfp(session; tol = 6, params...)
    nothing
end

@profview proftest()
