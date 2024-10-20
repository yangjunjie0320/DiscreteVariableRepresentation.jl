abstract type AbstractBasis end
abstract type AbstractPotential end

struct DiscreteVariableRepresentationProblem{T<:Real, B<:AbstractBasis}
    xg::Vector{T}
    vg::Vector{T}
    dx::T
    mu::T

    function DiscreteVariableRepresentationProblem(xg::Vector{T}, mu::T, b::B, v::AbstractPotential) where {T, B}
        ng = length(xg)
        @assert ng > 1
        @assert size(xg) == (ng,)
        
        dx = xg[2] - xg[1]
        @assert dx > 0.0

        for i in 2:ng
            xi = xg[i]; xj = xg[i-1]
            @assert isapprox(xi - xj, dx)
        end

        vg = [potential(v, x) for x in xg]
        @assert size(vg) == (ng,)
        
        new{T, B}(xg, vg, dx, mu)
    end
end

const DVR = DiscreteVariableRepresentationProblem

function hamiltonian(dvr_obj::DVR)
    hm = diagm(dvr_obj.vg)
    hm += kinetic(dvr_obj)

    ng = length(dvr_obj.xg)
    @assert size(hm) == (ng, ng)
    return hm
end

function solve(dvr_obj::DVR, nroots::Int=5)
    hm = hamiltonian(dvr_obj)
    res = eigen(Hermitian(hm), 1:nroots)
    e = res.values
    u = res.vectors
    return e, u
end

function plot(dvr_obj::DVR{T, B}, e::Vector{T}, u::Matrix{T}; nroots::Int=5, 
    xmin::Union{T, Nothing}=nothing, xmax::Union{T, Nothing}=nothing,
    ymin::Union{T, Nothing}=nothing, ymax::Union{T, Nothing}=nothing
    ) where {T, B}

    xg = dvr_obj.xg
    vg = dvr_obj.vg

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

    Plots.plot!(p, xlims=(xmin, xmax), ylims=(ymin, ymax))    
    return p
end
