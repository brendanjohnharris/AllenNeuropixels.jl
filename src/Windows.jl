function _slidingwindow(x::AbstractVector, width::Int=length(x)÷20, overlap::Int=0; tail=false)
    buf = 0:(width-overlap):(length(x)-width)
    𝓍 = [x[1+a:width+a] for a in buf]
    if tail == true
        thetail = x[last(buf)+width+1:end]
        length(thetail) > 0 && push!(𝓍, thetail)
    elseif tail == :overlap
        thetail = x[end-width+1:end]
        length(thetail) > 0 && push!(𝓍, thetail)
    end
    return 𝓍
end
function partition(x::AbstractVector, width::Int=50000)
    𝓍 = _slidingwindow(x, width, 0, tail=true)
end
function slidingwindow(x::AbstractVector, windowFunc::Function=rect, kwargs...)
    # windowFunc must be vectorised
    # Construct sub-timeseries
    X = hcat(_slidingwindow(x; kwargs...)...)
    idxs = 1:size(X, 1):(size(X, 2)*size(X, 1)+1)
    # Apply the window function to each sub-series
    mapslices(windowFunc, X, dims=1)
end
function slidingwindow(width::Int, args...)
    f(x) = slidingwindow(x, width, args...)
end
export slidingwindow

rect = identity
