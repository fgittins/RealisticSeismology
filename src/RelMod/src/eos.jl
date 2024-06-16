"""
Polytropic equations of state.

# Structs
- `Polytrope`
- `EnergyPolytrope`
- `StratifiedEnergyPolytrope`

# Functions
- `ε`
- `Γ`
- `Γ₁`

# Notes
Assumes geometric units, where G = c = 1.
"""

abstract type Barotrope end

"""
Adiabatic index [dimensionless] as function of pressure [km^-2].
"""
Γ₁(eos::Barotrope, p) = Γ(eos, p)

"""
Polytropic equation of state in terms of rest-mass density.

# Fields
- `n`: poltropic index.
- `K`: polytropic constant [km^(2/n)].
"""
struct Polytrope{T₁ <: Real, T₂ <: Real} <: Barotrope
    n::T₁
    K::T₂
end

"Energy density [km^-2] as function of pressure [km^-2]."
ε(eos::Polytrope, p) = (p / eos.K)^(eos.n/(eos.n + 1)) + eos.n*p

"""Background index [dimensionless] as function of pressure [km^-2]."""
Γ(eos::Polytrope, p) = 1 + 1/eos.n

"""
Polytropic equation of state in terms of energy density.

# Fields
- `n`: poltropic index.
- `K`: polytropic constant [km^(2/n)].
"""
struct EnergyPolytrope{T₁ <: Real, T₂ <: Real} <: Barotrope
    n::T₁
    K::T₂
end

ε(eos::EnergyPolytrope, p) = (p / eos.K)^(eos.n/(eos.n + 1))
Γ(eos::EnergyPolytrope, p) = (1 + 1/eos.n)*(1 + eos.K*(p / eos.K)^(1/(eos.n + 1)))

"""
Polytropic equation of state in terms of energy density with stratification.

# Fields
- `n`: poltropic index.
- `K`: polytropic constant [km^(2/n)].
- `frac`: ratio of adiabatic index to background index [dimensionless].
"""
struct StratifiedEnergyPolytrope{T₁ <: Real, T₂ <: Real, T₃ <: Real}
    n::T₁
    K::T₂
    frac::T₃
end

ε(eos::StratifiedEnergyPolytrope, p) = (p / eos.K)^(eos.n/(eos.n + 1))
Γ(eos::StratifiedEnergyPolytrope, p) = (1 + 1/eos.n)*(1 + eos.K*(p / eos.K)^(1/(eos.n + 1)))
Γ₁(eos::StratifiedEnergyPolytrope, p) = eos.frac*Γ(eos, p)
