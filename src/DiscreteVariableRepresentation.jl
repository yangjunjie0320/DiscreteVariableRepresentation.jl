module DiscreteVariableRepresentation

    using LinearAlgebra
    using Plots, Printf

    export SincBasis, DVR
    export solve, plot
    export SimpleHarmonicOscillator
    export MorseOscillator
    export SquareWell
    export potential

    include("potential.jl")
    include("dvr.jl")
    include("sinc.jl")
end
