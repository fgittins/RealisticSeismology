"""
Defines functions to calculate star with constant redshifted temperature.

# Structs
- `RedStar`

# Functions
- `structure`
- `integrate_star`
- `shooting`
- `RedStar`
- `Base.show`
"""

"""
Relativistic equations of stellar structure with constant redshifted
temperature.
"""
function structure(y, param, r)
    m, p, ν = y
    eos, pf, Tᵣ = param

    if p < 0
        εᵣ = zero(p)
    else
        εᵣ = ε(eos, Tᵣ/exp(ν/2), p)
    end

    dm_dr = 4*π*r^2*εᵣ
    dν_dr = 2*(m + 4*π*r^3*p)/(r*(r - 2*m))
    dp_dr = - (εᵣ + p)*dν_dr/2
    SA[dm_dr, dp_dr, dν_dr]
end

"Integrate stellar-structure equations with constant redshifted temperature."
function integrate_star(eos, Tᵣ, pc, pf, νc; save_everystep=false)
    r₀ = 1e-5
    εc = ε(eos, Tᵣ/exp(νc/2), pc)

    p₂ = - 4*π/3*(εc + pc)*(εc + 3*pc)
    ν₂ = 8*π/3*(εc + pc)
    m₃ = 4*π*εc

    m₀ = 0 + 1/3*r₀^3*m₃
    p₀ = pc + 1/2*r₀^2*p₂
    ν₀ = νc + 1/2*r₀^2*ν₂

    prob = ODEProblem{false}(structure, SA[m₀, p₀, ν₀], (r₀, 100),
                             (eos, pf, Tᵣ))
    sol = solve(prob, Vern9();
                callback=RelMod.cb, abstol=1e-10, reltol=1e-10,
                save_everystep=save_everystep)
end

"Shoot to determine central metric potential `νc` [dimensionless]."
function shooting(νc, param)
    eos, Tᵣ, pc = param

    pf = exp(eos.itp₁(log(Tᵣ/exp(νc/2)), eos.ymin))/pressure_geometric_to_natural

    sol = integrate_star(eos, Tᵣ, pc, pf, νc)

    R = sol.t[end]
    M = sol[1, end]

    νsol = @view sol[3, :]

    νsol[end] - log(1 - 2*M/R)
end

"""
Relativistic star with constant redshifted temperature.

# Fields
- `eos`: equation of state.
- `pc`: central pressure [km^-2].
- `νc`: central metric potential [dimensionless].
- `εc`: central energy density [km^-2].
- `r₀`: smallest radius [km].
- `R`: stellar radius [km].
- `M`: stellar mass [km].
- `sol`: numerical solution.
- `Tᵣ`: redshifted temperature [MeV].
"""
struct RedStar{EOSType,
               T₁ <: Real, T₂ <: Real, T₃ <: Real, T₄ <: Real, T₅ <: Real,
               T₆ <: Real,
               SolType,
               T₇ <: Real}
    eos::EOSType
    pc::T₁
    νc::T₂
    εc::T₃
    r₀::T₄
    R::T₅
    M::T₆
    sol::SolType
    Tᵣ::T₇          # Redshifted temperature [MeV]
end

"""
Relativistic star with constant redshifted temperature.

# Arguments
- `eos`: equation of state.
- `Tᵣ`: redshifted temperature [MeV].
- `pc`: central pressure [km^-2].
"""
function RedStar(eos, Tᵣ, pc)
    ν₀min, ν₀max = -1.0, 0.0
    prob = IntervalNonlinearProblem{false}(shooting, (ν₀min, ν₀max),
                                           (eos, Tᵣ, pc))
    sol = solve(prob, ITP())
    νc = sol.u

    pf = exp(eos.itp₁(log(Tᵣ/exp(νc/2)), eos.ymin))/pressure_geometric_to_natural
    sol = integrate_star(eos, Tᵣ, pc, pf, νc; save_everystep=true)

    r₀ = sol.t[1]
    ν₀ = sol[3, 1]
    R = sol.t[end]
    M = sol[1, end]

    εc = ε(eos, Tᵣ/exp(νc/2), pc)

    ν₂ = 8*π/3*(εc + pc)
    νc = ν₀ - 1/2*r₀^2*ν₂

    star = RedStar(eos, pc, νc, εc, r₀, R, M, sol, Tᵣ)
    star
end

function Base.show(io::IO, ::MIME"text/plain", star::RedStar)
    println(io, "Star with M = $(star.M) km, R = $(star.R) km")
    println(io, "    pc = $(star.pc) km^-2")
    println(io, "    νc = $(star.νc)")
    print(io, "    εc = $(star.εc) km^-2")
end
