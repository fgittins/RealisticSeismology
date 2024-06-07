"""
    RelMod

Calculates oscillation modes of relativistic star.
"""
module RelMod

export Polytrope, EnergyPolytrope, StratifiedEnergyPolytrope
export Star
export solve_eigenfrequency, Muller, Simplex

using Optim
using OrdinaryDiffEq
using Roots
using StaticArrays

include("eos.jl")
include("mode.jl")
include("perturbations_exterior.jl")
include("perturbations_interior.jl")
include("star.jl")

end # module RelMod
