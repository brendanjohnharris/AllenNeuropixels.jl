using Distances
import Clustering: hclust
using Random
using SparseArrays

function lfpspikesimilarity(X::LFPVector, Y::Vector; normalize = identity)
    ts = Interval(extrema(dims(X, 𝑡))...)
    Y = Y[Y .∈ (ts,)]
    X = normalize(X)
    return mean(X[𝑡(Near(Y))]) # Mean value of the LFP at each spike time
end

function lfpspikesimilarity(session, probeid, X::LFPMatrix, Y::Dict; normalize = identity,
                            kwargs...)
    X = normalize(X)
    spikes = values(Y) |> collect
    units = keys(Y) |> collect
    # * First, match up the units to their nearest LFPs
    depths = getchanneldepths(session, probeid, X)
    idxs = sortperm(depths)
    X = X[:, idxs]
    depths = depths[idxs]
    X = ToolsArray(X, (dims(X, 𝑡), Depth(depths)))
    unitdepths = getunitdepths(session, probeid, units)
    X = X[Depth(Near(unitdepths))] # Sorted in order of spikes
    sims = [lfpspikesimilarity(X[:, i], spikes[i]; kwargs...) for i in eachindex(spikes)]
    return Dict(units .=> sims)
end

function pairspikelfp(session, structure, spikes, X::LFPMatrix)
    probeid = getprobe(session, structure)
    units = keys(spikes) |> collect
    spikes = values(spikes) |> collect
    idxs = AllenNeuropixelsBase.isvalidunit(session, units)
    units = units[idxs]
    spikes = spikes[idxs]
    # * First, match up the units to their nearest LFPs
    depths = getchanneldepths(session, probeid, X)
    idxs = sortperm(depths)
    X = X[:, idxs]
    depths = depths[idxs]
    X = ToolsArray(X, (dims(X, 𝑡), Depth(depths)))
    unitdepths = getunitdepths(session, probeid, units)
    X = X[Depth(Near(unitdepths))] # Sorted in order of spikes
    int = Interval(extrema(dims(X, 𝑡))...)
    spikes = [s[s .∈ (int,)] for s in spikes]
    idxs = .!isempty.(spikes)
    spikes = spikes[idxs]
    units = units[idxs]
    X = X[:, idxs]
    return (units, spikes, X)
end

function pairspikelfp(session, structure::String, X::LFPMatrix)
    spikes = getspiketimes(session, structure)
    pairspikelfp(spikes, X)
end

function pairspikelfp(params::NamedTuple, args...)
    session = Session(params[:sessionid])
    structure = params[:structure]
    return pairspikelfp(session, structure, args...)
end

mape(x, y) = sum(1.0 .- abs.((x .- y) ./ x)) / length(x)
sse(x, y) = sum((x .- y) .^ 2)
msse(x, y) = std((x .- y) .^ 2)
R²(x, y) = 1.0 - sse(x, y) / sse(x, mean(x))
absoluteR²(x, y) = R²(x, y) > 0 ? R²(x, y) : R²(x, -y) # If the predictions are negatively correlated to x

function predictionerror(x, y, M::CCA; metric = R²)

    # * Forward prediction, x -> y
    zx = predict(M, x, :x)
    zy = zx # ? The assumption that the low-dim shared space maximises correlation between the two latent variables
    ŷ = (M.yproj)' \ zy .+ M.ymean

    # * Reverse prediction, y -> x
    zy = predict(M, y, :y)
    zx = zy
    x̂ = (M.xproj)' \ zx .+ M.xmean
    return metric(x, x̂), metric(y, ŷ)
end

"""
Calculate the prediction error of variables 1 to variables 2 and vice versa.
The output contains prediction errors.
If used with spike matrices probably want to transpose those
"""
function predictionerror(x, y; metric = R², model = MultivariateStats.CCA,
                         maxdim = min(size(x, 1), size(y, 1)), kwargs...)
    x = collect(x)
    y = collect(y)
    N = size(x, 2)
    nₜ = N ÷ 5
    iₜ = fill(false, N)
    iₜ[randperm(N)[1:nₜ]] .= true
    xₜ = x[:, iₜ]
    yₜ = y[:, iₜ]
    x = x[:, .!iₜ]
    y = y[:, .!iₜ]

    if model == MultivariateStats.CCA
        Δ = [predictionerror(x, y, fit(model, x, y; outdim = n, kwargs...); metric)
             for n in 1:maxdim]
    else
        Δ = [(metric(xₜ, predict(fit(model, y, x; r = n, kwargs...), yₜ)), # bckwd
              metric(yₜ, predict(fit(model, x, y; r = n, kwargs...), xₜ)))
             for n in 1:maxdim]
    end
    return first.(Δ), last.(Δ)
end
