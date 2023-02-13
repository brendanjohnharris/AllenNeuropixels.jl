function getrunningspeed(S::AbstractSession)
    f = h5open(getfile(S), "r")
    r = f["processing"]["running"]["running_speed"]["data"] |> read
    ts = f["processing"]["running"]["running_speed"]["timestamps"] |> read
    return DimArray(r, (Ti(ts),); metadata=Dict(:sessionid=>getid(S)))
end

function smoothrunningspeed(S::AbstractSession; windowfunc=hanning, window=1)
    r = getrunningspeed(S)
    ts = r.dims[1]
    dt = mean(diff(collect(r.dims[1])))
    n = ceil(Int, window/dt/2)*2
    x = DSP.conv(r, windowfunc(n))
    x = x[n÷2:end-n÷2]
    return DimArray(x, (Ti(collect(ts)),), metadata=r.metadata)
end

smoothrunningspeed(s::Integer; kwargs...) = smoothrunningspeed(Session(s); kwargs...)
