"""
x and y are vectors containing spike times
"""
function rbfdistance(x::AbstractVector, y::AbstractVector; σ, δ=3σ)
    # We will make the assumption that the scale (and amplitude) parameter is the same for every spike. I.e., all spikes are identical.
    # integral(exp(-1/2 ((t - t_1)/s)^2) exp(-1/2 ((t - t_2)/s)^2))/((s sqrt(2 π)) (s sqrt(2 π))) dt between - infinity and infinity
    # See ./test/spikeclusterintegral.jl
    # φ(tᵢ; σ=1) = t -> exp(-((t - tᵢ)/σ)^2/2)/(σ*sqrt(2π))
    # ! This is UNNORMALIZED. Slightly valid if the spikes are further apart than the kernel width.
    t⃗ = Iterators.product(x, y) |> collect |> vec # Pairs of spikes
    filter!(x->abs(reduce(-, x))<δ, t⃗)
    f = ((t₁, t₂),) -> exp(-(t₁ - t₂)^2)/(σ*sqrt(2π))
    I = f.(t⃗) |> sum
end

function rbfdistance(x::AbstractSparseDimArray, y::AbstractSparseDimArray; kwargs...)
    ix, iy = (x, y) .|> SparseVector .|> findnz .|> first
    x = dims(x, Ti)[ix] |> collect
    y = dims(y, Ti)[iy] |> collect
    rbfdistance(x, y; kwargs...)
end
