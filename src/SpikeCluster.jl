"""
x and y are vectors containing spike times
"""
function rbfdistance(x::Vector, y::Vector; σ, δ=3σ)
    φ(tᵢ; σ=1) = t -> exp(-((t - tᵢ)/σ)^2/2)/(σ*sqrt(2π))
    t⃗ = Iterators.product(x, y) # Pairs of neurons
    filter!(x->abs(diff(x))>δ, t⃗)
end
