module DiscreteVariableRepresentation

    using LinearAlgebra
    using Plots, Printf

    export Sinc, DVR
    export solve, plot
    export SimpleHarmonicOscillator
    export Morse, SquareWell
    export DoubleWell

    include("dvr.jl")
    include("potential.jl")
    include("basis.jl")
end
