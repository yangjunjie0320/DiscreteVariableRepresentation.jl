abstract type AbstractPotential end

struct SimpleHarmonicOscillator <: AbstractPotential 
    k::Real
    x0::Real
end

function potential(v::SimpleHarmonicOscillator, x::Real)
    k = v.k
    x0 = v.x0
    return k * (x - x0)^2 / 2
end

struct MorseOscillator <: AbstractPotential
    de::Real
    re::Real
    beta::Real
end

function potential(v::MorseOscillator, x::Real)
    de = v.de
    re = v.re
    beta = v.beta
    return de * (1 - exp(-beta * (x - re)))^2
end

struct SquareWell <: AbstractPotential
    v::Real
    x0::Real
    l::Real
end

function potential(v::SquareWell, x::Real)
    x0 = v.x0
    l = v.l
    return abs(x - x0) > l/2 ? v.v : 0.0
end
