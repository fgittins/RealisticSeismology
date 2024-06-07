"""
Generating stellar models for given baryon masses.

# Functions
- `baryon_mass`
- `calculate_baryon_mass`
- `f`
- `calculate_star`

# Notes
Extends `RelMod` package.
"""

using OrdinaryDiffEq
using SimpleNonlinearSolve

function baryon_mass(y, param, r)
    mb = y
    star = param

    mᵣ = RelMod.m(star, r)
    pᵣ = RelMod.p(star, r)

    nbᵣ = nb(star.eos, pressure_geometric_to_natural*pᵣ)
    ρ = nbᵣ*1e54*mₙ/mass_geometric_to_CGS

    dmb_dr = 4*π*r^2*ρ/√(1 - 2*mᵣ/r)
    dmb_dr
end

function calculate_baryon_mass(star)
    nbc = nb(star.eos, pressure_geometric_to_natural*star.pc)
    ρc = nbc*1e54*mₙ/mass_geometric_to_CGS

    mb₃ = 4*π*ρc

    mb₀ = 0 + 1/3*star.r₀^3*mb₃

    prob = ODEProblem{false}(baryon_mass, mb₀, (star.r₀, star.R), star)
    sol = solve(prob, Tsit5();
                abstol=1e-10, reltol=1e-10, save_everystep=false)

    Mb = sol[end]
    Mb
end

function f(pc, param)
    eos, Mbtarget, pf = param

    star = Star(eos, pc, pf)
    Mb = calculate_baryon_mass(star)
    Mb - Mbtarget
end

function calculate_star(eos, Mb, pf)
    prob = IntervalNonlinearProblem{false}(f, (1e-3, 1e-5), (eos, Mb, pf))
    sol = solve(prob, ITP())
    pc = sol.u

    star = Star(eos, pc, pf)
    star
end

calculate_star(eos, Mb) = calculate_star(eos, Mb, 0)
