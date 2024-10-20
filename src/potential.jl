abstract type AbstractPotential end

struct SquareWell <: AbstractPotential
    x1::Real
    x2::Real
    depth::Real

    function SquareWell(;x1::Real=-0.5, x2::Real=0.5, depth::Real=1.0)
        @assert x1 < x2
        @assert depth > 0
        new(x1, x2, depth)
    end
end

function poten(dvr_obj::DVR{T, B, SquareWell}) where {T, B}
    d = dvr_obj.v.depth::T
    x1 = dvr_obj.v.x1::T
    x2 = dvr_obj.v.x2::T

    xx = dvr_obj.b.xg
    yy = ones(T, length(xx)) * d
    yy[(xx .> x1) .& (xx .< x2)] .= 0.0
    return diagm(yy::Vector{T})
end

struct DoubleWell <: AbstractPotential
    x1::Real
    x2::Real
    x3::Real
    x4::Real
    depth::Real

    function DoubleWell(;x1::Real=-2.0, x2::Real=-1.0, x3::Real=1.0, x4::Real=2.0, depth::Real=1.0)
        @assert x1 < x2 < x3 < x4
        @assert depth > 0
        new(x1, x2, x3, x4, depth)
    end
end

function poten(dvr_obj::DVR{T, B, DoubleWell}) where {T, B}
    d  = dvr_obj.v.depth::T
    x1 = dvr_obj.v.x1::T
    x2 = dvr_obj.v.x2::T
    x3 = dvr_obj.v.x3::T
    x4 = dvr_obj.v.x4::T

    xx = dvr_obj.b.xg
    yy = ones(T, length(xx)) * d
    yy[(xx .> x1) .& (xx .< x2)] .= 0.0
    yy[(xx .> x3) .& (xx .< x4)] .= 0.0
    return diagm(yy::Vector{T})
end

struct SimpleHarmonicOscillator <: AbstractPotential 
    k::Real
    x0::Real

    function SimpleHarmonicOscillator(;k::Real=1.0, x0::Real=0.0)
        @assert k > 0
        new(k, x0)
    end
end

function poten(dvr_obj::DVR{T, B, SimpleHarmonicOscillator}) where {T, B}
    k = dvr_obj.v.k::T
    x0 = dvr_obj.v.x0::T

    xx = dvr_obj.b.xg .- x0
    yy = k * xx .^ 2 / 2
    return diagm(yy::Vector{T})
end

struct Morse <: AbstractPotential
    d::Real
    a::Real
    x0::Real

    function Morse(;d::Real=1.0, a::Real=1.0, x0::Real=0.0)
        @assert d > 0
        @assert a > 0
        new(d, a, x0)
    end
end

function poten(dvr_obj::DVR{T, B, Morse}) where {T, B}
    d = dvr_obj.v.d::T
    a = dvr_obj.v.a::T
    x0 = dvr_obj.v.x0::T

    xx = dvr_obj.b.xg .- x0
    yy = d * (1 .- exp.(-a * xx)) .^ 2 .- d
    return diagm(yy::Vector{T})
end
