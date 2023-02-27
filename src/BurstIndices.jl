function convolvespikes(dt::Number, s::Vector; window, windowfunc=hanning)
    n = ceil(Int, window/dt/2)*2 + 1 # Window width, always odd
    w = windowfunc(n)
    w = w./sum(w)./dt # Normalize the window so that each spike contributes '1' rate. Is this right?
    # Assume the window is dt*n wide
    ts = (minimum(s) - dt*(n-1)/2):dt:(maximum(s) + dt*(n-1)/2)
    x = DimArray(zeros(length(ts)), (Ti(ts),))
    for t in s
        t = ts[findmin(abs.(ts .- t))[2]] # Centre of the window
        windowinds = (t-dt*(n-1)//2):dt:(t+dt*(n-1)//2 + dt/10)
        x[Ti(Near(windowinds))] .+= w
    end
    return x
end

function burstratetrace(B::BurstVector; window=5, kwargs...) # Default 5s window
    δ = mean(dt.(B)) # timestep for resulting trace
    s = peaktime.(B)
    s = convolvespikes(δ, s; window, kwargs...)
end
