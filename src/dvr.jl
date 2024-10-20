abstract type AbstractBasis end
abstract type AbstractPotential end

struct DiscreteVariableRepresentationProblem{T<:Real, B<:AbstractBasis, V<:AbstractPotential}
    v::V # potential
    b::B # basis
    mu::T

    function DiscreteVariableRepresentationProblem(v::V, b::B, mu::T=1.0) where {T, B, V}
        @assert mu > 0.0
        new{T, B, V}(v, b, mu)
    end
end

function DiscreteVariableRepresentationProblem(xg::Vector{T}; v::Union{V, Nothing}=nothing, b::Union{B, Nothing}=nothing, mu::T=1.0) where {T, V<:AbstractPotential, B<:AbstractBasis}
    b = b === nothing ? Sinc(xg) : b
    v = v === nothing ? SimpleHarmonicOscillator() : v
    return DiscreteVariableRepresentationProblem(v, b, mu)
end

const DVR = DiscreteVariableRepresentationProblem

function hamiltonian(dvr_obj::DVR{T, B, V}) where {T, B, V}
    vm = poten(dvr_obj)
    km = kinet(dvr_obj)
    hm = vm + km

    ng = dvr_obj.b.ng
    @assert size(hm) == (ng, ng)

    @assert ishermitian(hm)
    return Hermitian(hm::Matrix{T})
end

function solve(dvr_obj::DVR, nroots::Int=5)
    nroots = min(nroots, dvr_obj.b.ng)::Int
    @assert nroots > 0

    hm = hamiltonian(dvr_obj)
    res = eigen(hm, 1:nroots)
    e = res.values
    u = res.vectors
    return e, u
end

function plot(
    dvr_obj::DVR{T, V, B}, e::Vector{T}, u::Matrix{T}; nroots::Int=5, 
    xmin::Union{T, Nothing}=nothing, xmax::Union{T, Nothing}=nothing,
    ymin::Union{T, Nothing}=nothing, ymax::Union{T, Nothing}=nothing
    ) where {T, V, B}

    xg = dvr_obj.b.xg
    vg = diag(poten(dvr_obj))

    nroots = min(nroots, size(u, 2))
    ee = e[1:nroots]
    uu = u[:, 1:nroots]

    xmin = xmin === nothing ? minimum(xg) : xmin
    xmax = xmax === nothing ? maximum(xg) : xmax
    ymin = ymin === nothing ? ceil(minimum(vg) - 1.0) : ymin
    ymax = ymax === nothing ? floor(max(maximum(abs.(uu) .+ ee'), maximum(vg)) + 1.0) : ymax

    p = Plots.plot(xg, vg)
    for i in 1:nroots
        ix = argmax(abs.(uu[:, i]))
        u0 = uu[:, i] / uu[ix, i] * abs(uu[ix, i])
        Plots.plot!(p, xg, u0 .+ ee[i])
    end

    Plots.plot!(p, xlims=(xmin, xmax), ylims=(ymin, ymax), legend=false)
    return p
end
