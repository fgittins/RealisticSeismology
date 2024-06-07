"""
    Redshift

Calculates oscillation modes of relativistic star at constant redshifted
temperature.
"""
module Redshift

export APR
export RedStar
export calculate_baryon_mass, calculate_star
export calculate_love_number
export solve_eigenfrequency, Muller, Simplex

using Interpolations
using OrdinaryDiffEq
using RelMod
using SimpleNonlinearSolve
using StaticArrays

import RelMod: integrate_interior, integrate_interior_low_frequency
import RelMod: solve_eigenfrequency, Muller, Simplex
import RelMod: spectrum

include("../../read_table.jl")
include("../../units.jl")

include("apr.jl")
include("star.jl")
include("baryon_mass.jl")
include("love.jl")
include("perturbations_interior.jl")

end # module Redshift
