"""
Structure of relativistic star.

# Structs
- `Star`

# Functions
- `structure!`
- `surface`
- `integrate_star`
- `Star`
- `m`
- `p`
- `ν`
- `Base.show`

# Constants
- `cb`

# Notes
Assumes geometric units, where G = c = 1.
"""

"Relativistic equations of stellar structure."
function structure!(dy_dr, y, param, r)
    m, p, ν = y
    eos, pf = param

    if p < 0
        εᵣ = zero(p)
    else
        εᵣ = ε(eos, p)
    end

    dy_dr[1] = dm_dr = 4*π*r^2*εᵣ
    dy_dr[3] = dν_dr = 2*(m + 4*π*r^3*p)/(r*(r - 2*m))
    dy_dr[2] = dp_dr = - (εᵣ + p)*dν_dr/2
    dy_dr
end

"Definition of stellar surface."
surface(y, r, integrator) = y[2] - integrator.p[2]

affect!(integrator) = terminate!(integrator)

const cb = ContinuousCallback(surface, affect!)

"Integrate stellar-structure equations."
function integrate_star(eos, pc, pf)
    r₀ = 1e-5
    εc = ε(eos, pc)

    p₂ = - 4*π/3*(εc + pc)*(εc + 3*pc)
    m₃ = 4*π*εc

    m₀ = 0 + 1/3*r₀^3*m₃
    p₀ = pc + 1/2*r₀^2*p₂

    prob = ODEProblem{true}(structure!, [m₀, p₀, 0], (r₀, 100), (eos, pf))
    sol = solve(prob, Vern9(); callback=cb, abstol=1e-10, reltol=1e-10)

    νsol = @view sol[3, :]
    R = sol.t[end]
    M = sol[1, end]

    νsol .= νsol .+ (log(1 - 2*M/R) - νsol[end])

    sol
end

"""
Relativistic star.

# Fields
- `eos`: equation of state.
- `pc`: central pressure [km^-2].
- `νc`: central metric potential [dimensionless].
- `εc`: central energy density [km^-2].
- `r₀`: smallest radius [km].
- `R`: stellar radius [km].
- `M`: stellar mass [km].
- `sol`: numerical solution.
"""
struct Star{EOSType,
            T₁ <: Real, T₂ <: Real, T₃ <: Real, T₄ <: Real, T₅ <: Real,
            T₆ <: Real,
            SolType}
    eos::EOSType
    pc::T₁
    νc::T₂
    εc::T₃
    r₀::T₄
    R::T₅
    M::T₆
    sol::SolType
end

"""
Relativistic star.

# Arguments
- `eos`: equation of state.
- `pc`: central pressure [km^-2].
- `pf`: surface pressure [km^-2].
"""
function Star(eos, pc, pf)
    sol = integrate_star(eos, pc, pf)

    r₀ = sol.t[1]
    ν₀ = sol[3, 1]
    R = sol.t[end]
    M = sol[1, end]

    εc = ε(eos, pc)

    ν₂ = 8*π/3*(εc + pc)
    νc = ν₀ - 1/2*r₀^2*ν₂

    Star(eos, pc, νc, εc, r₀, R, M, sol)
end

Star(eos, pc) = Star(eos, pc, 0)

"Mass [km] as function of radius [km]."
m(star, r) = star.sol(r; idxs=1)

"Pressure [km^-2] as function of radius [km]."
p(star, r) = star.sol(r; idxs=2)

"Metric potential [dimensionless] as function of radius [km]."
ν(star, r) = star.sol(r; idxs=3)

function Base.show(io::IO, ::MIME"text/plain", star::Star)
    println(io, "Star with M = $(star.M) km, R = $(star.R) km")
    println(io, "    pc = $(star.pc) km^-2")
    println(io, "    νc = $(star.νc)")
    print(io, "    εc = $(star.εc) km^-2")
end
