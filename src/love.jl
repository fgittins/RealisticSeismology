"""
Calculating quadrupolar tidal Love number using formulation in Ref. [1].

# Functions
- `static_perturbations_interior`
- `calculate_love_number`

# Notes
Extends `RelMod` package.

# References
[1] Hinderer, "Tidal Love Numbers of Neutron Stars," Astrophys. J. 677 (2),
    1216 (2008).
"""

using OrdinaryDiffEq

function static_perturbations_interior(y, param, r)
    background, l = param

    mᵣ = RelMod.m(background, r)
    pᵣ = RelMod.p(background, r)

    εᵣ = RelMod.ε(background.eos, pᵣ)
    Γᵣ = RelMod.Γ(background.eos, pᵣ)

    expλ = 1/(1 - 2*mᵣ/r)
    dν_dr = 2*(mᵣ + 4*π*r^3*pᵣ)/(r*(r - 2*mᵣ))

    dy_dr = (- (y - 1)*y - (2 + expλ*(2*mᵣ/r + 4*π*r^2*(pᵣ - εᵣ)))*y
             + l*(l + 1)*expλ
             - 4*π*r^2*expλ*(5*εᵣ + 9*pᵣ + (εᵣ + pᵣ)^2/(Γᵣ*pᵣ))
             + (r*dν_dr)^2)/r
    dy_dr
end

function calculate_love_number(background)
    l = 2
    y₀ = l

    prob = ODEProblem{false}(static_perturbations_interior, Float64(y₀),
                             (background.r₀, background.R), (background, l))
    sol = solve(prob, Tsit5();
                abstol=1e-10, reltol=1e-10, save_everystep=false)

    C = background.M / background.R
    Y = sol[end]

    k₂ = (8*C^5/5*(1 - 2*C)^2*(2 + 2*C*(Y - 1) - Y)
          *(2*C*(6 - 3*Y + 3*C*(5*Y - 8))
            + 4*C^3*(13 - 11*Y + C*(3*Y - 2) + 2*C^2*(1 + Y))
            + 3*(1 - 2*C)^2*(2 - Y + 2*C*(Y - 1))*log(1 - 2*C))^(-1))
    k₂
end
