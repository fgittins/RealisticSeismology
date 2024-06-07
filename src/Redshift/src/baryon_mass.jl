"""
Generating stellar models for given baryon masses.

# Functions
- `baryon_mass`
- `calculate_baryon_mass`
- `f`
- `calculate_star`
"""

function baryon_mass(y, param, r)
    mb = y
    star = param

    mᵣ = RelMod.m(star, r)
    pᵣ = RelMod.p(star, r)
    νᵣ = RelMod.ν(star, r)

    nbᵣ = nb(star.eos, star.Tᵣ/exp(νᵣ/2), pressure_geometric_to_natural*pᵣ)
    ρ = nbᵣ*1e54*mₙ/mass_geometric_to_CGS

    dmb_dr = 4*π*r^2*ρ/√(1 - 2*mᵣ/r)
    dmb_dr
end

function calculate_baryon_mass(star)
    r₀ = star.sol.t[1]
    nbc = nb(star.eos, star.Tᵣ/exp(star.νc/2), pressure_geometric_to_natural*star.pc)
    ρc = nbc*1e54*mₙ/mass_geometric_to_CGS

    mb₃ = 4*π*ρc

    mb₀ = 0 + 1/3*r₀^3*mb₃

    prob = ODEProblem{false}(baryon_mass, mb₀, (r₀, star.R), star)
    sol = solve(prob, Tsit5();
                abstol=1e-10, reltol=1e-10, save_everystep=false)

    Mb = sol[end]
    Mb
end

function f(pc, param)
    eos, Tᵣ, Mbtarget = param

    star = RedStar(eos, Tᵣ, pc)
    Mb = calculate_baryon_mass(star)
    Mb - Mbtarget
end

function calculate_star(eos, Tᵣ, Mb)
    prob = IntervalNonlinearProblem{false}(f, (8e-4, 1e-5), (eos, Tᵣ, Mb))
    sol = solve(prob, ITP())
    pc = sol.u

    star = RedStar(eos, Tᵣ, pc)
    star
end
