using Dierckx

function getrunningspeed(S::AbstractSession)
    f = h5open(getfile(S), "r")
    r = f["processing"]["running"]["running_speed"]["data"] |> read
    ts = f["processing"]["running"]["running_speed"]["timestamps"] |> read
    return DimArray(r, (𝑡(ts),); metadata = Dict(:sessionid => getid(S)))
end

function smoothrunningspeed(S::AbstractSession; windowfunc = hanning, window = 1)
    r = getrunningspeed(S)
    ts = r.dims[1]
    dt = mean(diff(collect(r.dims[1])))
    n = ceil(Int, window / dt / 2) * 2
    w = windowfunc(n)
    w = w ./ sum(w) # Normalized wavelet. Outputs after convolution become averages
    x = DSP.conv(r, w)
    x = x[(n ÷ 2):(end - n ÷ 2)]
    return DimArray(x, (𝑡(collect(ts)),), metadata = r.metadata)
end

smoothrunningspeed(s::Integer; kwargs...) = smoothrunningspeed(Session(s); kwargs...)

function stimulustrace(S::AbstractSession, feature, times)
    times = Interval(extrema(times)...)
    df = getstimuli(S, times)
    s = df[:, feature]
    x = zeros(length(s))
    x[s .== "null"] .= NaN
    idxs = s .!= "null"
    x[idxs] .= Meta.parse.(s[idxs])
    t = mean.(zip(df.start_time, df.stop_time))
    return DimArray(x, (𝑡(t),))
end

function stimulustrace(S::AbstractSession, feature, times::AbstractRange)
    s = stimulustrace(S, feature, extrema(times))
    interpmatch(s, times)
end

function stimulustrace(S::AbstractSession, feature, x::LFPVector)
    t = dims(x, 𝑡)
    stimulustrace(S, feature, t)
end

function interpmatch(x::LFPVector, ts::AbstractRange)
    f = Spline1D(collect(dims(x, 𝑡)), x)
    x̂ = f(ts)
    return DimArray(x̂, (𝑡(ts),); metadata = x.metadata)
end

"""
Match the time indices of the first input DimVector to the second by interpolating the first
"""
function interpmatch(x::LFPVector, y::LFPVector)
    ts = dims(y, 𝑡).val.data
    interpmatch(x, ts)
end

# function interpcorr(x::LFPVector, y::LFPVector)
# 	x = interpmatch(x, y)
# 	r = corspearman(x, y)
# end
