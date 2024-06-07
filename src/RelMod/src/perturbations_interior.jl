"""
Define functions for mode perturbations in interior of relativistic star.

Based on equations derived in Ref. [1]. Low-frequency formulation comes from
Ref. [2].

# Functions
- `perturbations_interior`
- `taylor_coefficients`
- `integrate_interior`
- `perturbations_interior_low_frequency`
- `taylor_coefficients_low_frequency`
- `integrate_interior_low_frequency`
- `X`

# Notes
Assumes geometric units, where G = c = 1.

# References
[1] Detweiler and Lindblom, "On the nonradial pulsations of general
    relativistic stellar models," Astrophys. J. 292, 12 (1985).
[2] Krüger, "Seismology of adolescent general relativistic neutron stars," PhD
    thesis, University of Southampton (2015).
"""

"Interior polar perturbation equations for relativistic star."
function perturbations_interior(y, param, r)
    H₁, K, W, X = y
    background, l, ω² = param

    # background
    mᵣ = m(background, r)
    pᵣ = p(background, r)
    νᵣ = ν(background, r)
    εᵣ = ε(background.eos, pᵣ)
    Γ₁ᵣ = Γ₁(background.eos, pᵣ)

    dm_dr = 4*π*r^2*εᵣ
    dν_dr = 2*(mᵣ + 4*π*r^3*pᵣ)/(r*(r - 2*mᵣ))
    dp_dr = - (εᵣ + pᵣ)*dν_dr/2

    expν = exp(νᵣ)
    expλ = 1 / (1 - 2*mᵣ/r)
    dλ_dr = 2*expλ*(dm_dr/r - mᵣ/r^2)

    d²ν_dr² = 2/(r*(r - 2*mᵣ))*(dm_dr + 4*π*r^2*(3*pᵣ + r*dp_dr)
                                + (mᵣ + r*dm_dr - r)*dν_dr)

    # perturbations
    H₀ = (
        8*π*r^3/√expν*X
        - (l*(l + 1)/2*(mᵣ + 4*π*r^3*pᵣ)
           - ω²*r^3/(expλ*expν))*H₁
        + ((l + 2)*(l - 1)/2*r - ω²*r^3/expν
           - expλ/r*(mᵣ + 4*π*r^3*pᵣ)*(3*mᵣ - r + 4*π*r^3*pᵣ))*K
    ) / (3*mᵣ + (l + 2)*(l - 1)/2*r + 4*π*r^3*pᵣ)
    V = √expν/ω²*(
        1/(εᵣ + pᵣ)*X
        - dν_dr/(2*r)*√(expν/expλ)*W
        - √expν/2*H₀
    )

    dH₁_dr = (
        - (l + 1 + 2*mᵣ*expλ/r + 4*π*r^2*expλ*(pᵣ - εᵣ))/r*H₁
        + expλ/r*(
            H₀
            + K
            - 16*π*(εᵣ + pᵣ)*V
        )
    )
    dK_dr = (
        1/r*H₀
        + l*(l + 1)/(2*r)*H₁
        - ((l + 1)/r - dν_dr/2)*K
        - 8*π*(εᵣ + pᵣ)*√expλ/r*W
    )
    dW_dr = (
        - (l + 1)/r*W
        + r*√expλ*(
            1/(Γ₁ᵣ*pᵣ*√expν)*X
            - l*(l + 1)/r^2*V
            + 1/2*H₀
            + K
        )
    )
    dX_dr = (
        - l/r*X
        + (εᵣ + pᵣ)*√expν*(
            (1/r - dν_dr/2)/2*H₀
            + (r*ω²/expν + l*(l + 1)/(2*r))/2*H₁
            + (3*dν_dr/2 - 1/r)/2*K
            - l*(l + 1)*dν_dr/(2*r^2)*V
            - (4*π*(εᵣ + pᵣ)*√expλ
               + ω²*√expλ/expν
               + (dν_dr*(dλ_dr/2 + 2/r)
                  - d²ν_dr²)/(2*√expλ))/r*W
        )
    )
    SA[dH₁_dr, dK_dr, dW_dr, dX_dr]
end

"Calculate coefficients for Taylor expansion near centre."
function taylor_coefficients(Kc, Wc, param)
    background, l, ω² = param

    pc = background.pc
    εc = background.εc
    expνc = exp(background.νc)
    Γc = Γ(background.eos, pc)
    Γ₁c = Γ₁(background.eos, pc)

    # coefficients for background quantities
    p₂ = - 4*π/3*(εc + pc)*(εc + 3*pc)
    ε₂ = p₂*(εc + pc)/(Γc*pc)
    ν₂ = 8*π/3*(εc + 3*pc)

    p₄ = - (2*π/5*(εc + pc)*(ε₂ + 5*p₂)
            + 2*π/3*(ε₂ + p₂)*(εc + 3*pc)
            + 32*π^2/9*εc*(εc + pc)*(εc + 3*pc))
    ν₄ = (4*π/5*(ε₂ + 5*p₂)
          + 64*π^2/9*εc*(εc + 3*pc))

    # zeroth-order coefficients
    H₁c = (
        2*l*Kc
        + 16*π*(εc + pc)*Wc
    )/(l*(l + 1))
    Xc = (εc + pc)*√expνc*(
        1/2*Kc
        + (ν₂/2 - ω²/(l*expνc))*Wc
    )

    # solve for second-order coefficients
    Q₀ = 4/((l + 2)*(l - 1))*(
        8*π/√expνc*Xc
        - (8*π/3*εc + ω²/expνc)*Kc
        - (2*π/3*l*(l + 1)*(εc + 3*pc) - ω²/expνc)*H₁c
    )
    Q₁ = 2/(l*(l + 1))*(
        1/(Γ₁c*pc*√expνc)*Xc + 3/2*Kc
        + 4*π/3*(l + 1)*εc*Wc
    )

    b = SA[
        (1/4*ν₂/√expνc*Xc + 1/4*(ε₂ + p₂)*Kc + 1/4*(εc + pc)*Q₀
         + 1/2*ω²*(εc + pc)/expνc*Q₁
         - (p₄ - 4*π/3*εc*p₂
            + ω²/(2*l)*(ε₂ + p₂ - (εc + pc)*ν₂)/expνc)*Wc)
        ###
        4*π/3*(εc + 3*pc)*Kc + 1/2*Q₀ - 4*π*(ε₂ + p₂ + 8*π/3*εc*(εc + pc))*Wc
        ###
        (4*π*(1/3*(2*l + 3)*εc - pc)*H₁c + 8*π/l*(ε₂ + p₂)*Wc
         - 8*π*(εc + pc)*Q₁ + 1/2*Q₀)
        ###
        (1/2*(ε₂ + p₂ + 1/2*(εc + pc)*ν₂)*l/(εc + pc)*Xc
         + (εc + pc)*√expνc*(
                1/2*ν₂*Kc + 1/4*Q₀ + 1/2*ω²/expνc*H₁c
                - 1/4*l*(l + 1)*ν₂*Q₁ + (
                                1/2*(l + 1)*ν₄
                                - 2*π*(ε₂ + p₂)
                                - 16*π^2/3*εc*(εc + pc)
                                + 1/2*(ν₄ - 4*π/3*εc*ν₂)
                                + 1/2*ω²/expνc*(ν₂ - 8*π/3*εc)
                                )*Wc))
    ]

    A = SA[
        0 -1/4*(εc + pc) 1/2*(p₂ + (εc + pc)*ω²*(l + 3)/(l*(l + 1)*expνc)) 1/(2*√expνc)
        ###
        -1/4*l*(l + 1) 1/2*(l + 2) 4*π*(εc + pc) 0
        ###
        1/2*(l + 3) -1 -8*π*(εc + pc)*(l + 3)/(l*(l + 1)) 0
        ###
        -1/8*l*(l + 1)*(εc + pc)*√expνc 0 -(εc + pc)*√expνc*(1/4*(l + 2)*ν₂ - 2*π*(εc + pc) - 1/2*ω²/expνc) 1/2*(l + 2)
    ]

    x = A \ b
    d²H₁_dr²c, d²K_dr²c, d²W_dr²c, d²X_dr²c = x
    H₁c, d²H₁_dr²c, d²K_dr²c, d²W_dr²c, Xc, d²X_dr²c
end

"Integrate perturbation equations in interior."
function integrate_interior(background, l, ω²; save_everystep=false)
    R = background.R
    # starting point
    r₀ = 1e-3
    # matching point
    rmatch = R/2

    # from centre
    # Solution 1
    Kc, Wc = one(ω²), zero(ω²)
    H₁c, d²H₁_dr²c, d²K_dr²c, d²W_dr²c, Xc, d²X_dr²c = taylor_coefficients(
            Kc, Wc, (background, l, ω²))
    H₁₀ = H₁c + 1/2*r₀^2*d²H₁_dr²c
    K₀ = Kc + 1/2*r₀^2*d²K_dr²c
    W₀ = Wc + 1/2*r₀^2*d²W_dr²c
    X₀ = Xc + 1/2*r₀^2*d²X_dr²c

    prob1 = ODEProblem{false}(perturbations_interior, SA[H₁₀, K₀, W₀, X₀],
                              (r₀, rmatch), (background, l, ω²))
    sol1 = solve(prob1, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # Solution 2
    Kc, Wc = zero(ω²), one(ω²)
    H₁c, d²H₁_dr²c, d²K_dr²c, d²W_dr²c, Xc, d²X_dr²c = taylor_coefficients(
            Kc, Wc, (background, l, ω²))
    H₁₀ = H₁c + 1/2*r₀^2*d²H₁_dr²c
    K₀ = Kc + 1/2*r₀^2*d²K_dr²c
    W₀ = Wc + 1/2*r₀^2*d²W_dr²c
    X₀ = Xc + 1/2*r₀^2*d²X_dr²c

    prob2 = ODEProblem{false}(perturbations_interior, SA[H₁₀, K₀, W₀, X₀],
                              (r₀, rmatch), (background, l, ω²))
    sol2 = solve(prob2, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # from surface
    # Solution 3
    H₁f, Kf, Wf = one(ω²), zero(ω²), zero(ω²)
    Xf = zero(ω²)

    prob3 = ODEProblem{false}(perturbations_interior, SA[H₁f, Kf, Wf, Xf],
                              (R, rmatch), (background, l, ω²))
    sol3 = solve(prob3, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # Solution 4
    H₁f, Kf, Wf = zero(ω²), one(ω²), zero(ω²)
    Xf = zero(ω²)

    prob4 = ODEProblem{false}(perturbations_interior, SA[H₁f, Kf, Wf, Xf],
                              (R, rmatch), (background, l, ω²))
    sol4 = solve(prob4, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # Solution 5
    H₁f, Kf, Wf = zero(ω²), zero(ω²), one(ω²)
    Xf = zero(ω²)

    prob5 = ODEProblem{false}(perturbations_interior, SA[H₁f, Kf, Wf, Xf],
                              (R, rmatch), (background, l, ω²))
    sol5 = solve(prob5, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # solve for coefficients of general solution
    # with W normalised to 1 at surface
    A = SA[0 0 sol3[3, 1] sol4[3, 1] sol5[3, 1]
         sol1[1, end] sol2[1, end] -sol3[1, end] -sol4[1, end] -sol5[1, end]
         sol1[2, end] sol2[2, end] -sol3[2, end] -sol4[2, end] -sol5[2, end]
         sol1[3, end] sol2[3, end] -sol3[3, end] -sol4[3, end] -sol5[3, end]
         sol1[4, end] sol2[4, end] -sol3[4, end] -sol4[4, end] -sol5[4, end]]
    b = SA[1, 0, 0, 0, 0]
    x = A \ b

    x, sol1, sol2, sol3, sol4, sol5
end

"""
Interior polar perturbation equations for relativistic star with formulation
appropriate for low-frequency oscillations.
"""
function perturbations_interior_low_frequency(y, param, r)
    H₁, K, W, V = y
    background, l, ω² = param

    # background
    mᵣ = m(background, r)
    pᵣ = p(background, r)
    νᵣ = ν(background, r)
    εᵣ = ε(background.eos, pᵣ)
    Γᵣ = Γ(background.eos, pᵣ)
    Γ₁ᵣ = Γ₁(background.eos, pᵣ)

    expλ = 1 / (1 - 2*mᵣ/r)
    expν = exp(νᵣ)

    dλ_dr = (1 - expλ)/r + 8*π*r*expλ*εᵣ
    dν_dr = (expλ - 1)/r + 8*π*r*expλ*pᵣ
    dp_dr = - (εᵣ + pᵣ)*dν_dr/2

    n = (l + 2)*(l - 1)/2
    A = (1/Γᵣ - 1/Γ₁ᵣ)*dp_dr/pᵣ

    # perturbations
    H₀ = (
        r^2/expλ*(ω²*r/expν - (n + 1)*dν_dr/2)*H₁
        + (n*r - ω²*r^3/expν - r^2*dν_dr/(4*expλ)*(r*dν_dr - 2))*K
        + 4*π*r^2*(εᵣ + pᵣ)*(
            dν_dr/√expλ*W
            + 2*r*ω²/expν*V
        )
    ) / ((n + 1)*r - r/(2*expλ)*(r*dλ_dr + 2))

    dH₁_dr = (
        ((dλ_dr - dν_dr)/2 - (l + 1)/r)*H₁
        + expλ/r*(
            H₀
            + K
            - 16*π*(εᵣ + pᵣ)*V
        )
    )
    dK_dr = (
        1/r*H₀
        + (n + 1)/r*H₁
        + (dν_dr/2 - (l + 1)/r)*K
        - 8*π*(εᵣ + pᵣ)*√expλ/r*W
    )
    dW_dr = (
        - ((l + 1)/r + dp_dr/(Γ₁ᵣ*pᵣ))*W
        + r*√expλ*(
            (εᵣ + pᵣ)/(Γ₁ᵣ*pᵣ)*(
                ω²/expν*V
                + 1/2*H₀
            )
            - 2*(n + 1)/r^2*V
            + 1/2*H₀
            + K
        )
    )
    dV_dr = (
        (- A + dν_dr - l/r)*V
        - expν/(2*ω²)*A*(
            H₀
            + dν_dr/(r*√expλ)*W
        )
        + r*H₁
        - √expλ/r*W
    )
    SA[dH₁_dr, dK_dr, dW_dr, dV_dr]
end

"""
Calculate coefficients for Taylor expansion near centre with formulation
appropriate to low-frequency oscillations.
"""
function taylor_coefficients_low_frequency(Kc, Wc, param)
    background, l, ω² = param

    pc = background.pc
    εc = background.εc
    expνc = exp(background.νc)
    Γc = Γ(background.eos, pc)
    Γ₁c = Γ₁(background.eos, pc)

    # coefficients for background quantities
    λ₂ = 16*π/3*εc
    ν₂ = 8*π/3*(εc + 3*pc)
    p₂ = - ν₂/2*(εc + pc)
    ε₂ = p₂*(εc + pc)/(Γc*pc)

    # zeroth-order coefficients
    H₁c = (
        2*l*Kc
        + 16*π*(εc + pc)*Wc
    ) / (l*(l + 1))
    Vc = -1/l*Wc

    n = (l + 2)*(l - 1)/2

    # solve for second-order coefficients
    b = SA[
        (ν₂ - λ₂)/2*H₁c - λ₂*Kc + 8*π*(ε₂ + p₂ + λ₂*(εc + pc))*Vc
        ###
        - ν₂/2*Kc + 2*π*(2*(ε₂ + p₂) + λ₂*(εc + pc))*Wc
        ###
        (- 1/2*(3*Γ₁c*pc + pc + εc)*Kc + ((n + 1)/2*λ₂*Γ₁c*pc - l*p₂
         - ω²/expνc*(εc + pc))*Vc)
        ###
        (2*H₁c - expνc/ω²*(ε₂/(εc + pc) - p₂/(Γ₁c*pc))*Kc
         + (l/2*λ₂ + 2*ν₂
            - (2 - l*expνc*ν₂/ω²)*(ε₂/(εc + pc)
                                   - p₂/(Γ₁c*pc)))*Vc)
        ###
        (((n + 1)/2*ν₂ - ω²/expνc)*H₁c + (ω²/expνc - ν₂/2)*Kc
         - 8*π*ω²/expνc*(εc + pc)*Vc + 8*π*p₂*Wc)
    ]

    A = SA[
        -(l + 3)/2 1/2 0 -8*π*(εc + pc) 1/2
        ###
        (n + 1)/2 -(l + 3)/2 -4*π*(εc + pc) 0 1/2
        ###
        0 0 -(l + 3)/2*Γ₁c*pc -(n + 1)*Γ₁c*pc 0
        ###
        0 0 1 l + 2 0
        ###
        0 n/2 0 0 -n/2
    ]

    x = A \ b
    d²H₁_dr²c, d²K_dr²c, d²W_dr²c, d²Vdr²c, d²H₀dr²c = x
    H₁c, d²H₁_dr²c, d²K_dr²c, d²W_dr²c, Vc, d²Vdr²c
end

"""
Integrate perturbation equations in interior with formulation appropriate for
low-frequency oscillations.
"""
function integrate_interior_low_frequency(background, l, ω²;
                                          save_everystep=false)
    R = background.R
    # starting point
    r₀ = 1e-3
    # matching point
    rmatch = R/2

    # from centre
    # Solution 1
    Kc, Wc = one(ω²), zero(ω²)
    (H₁c, d²H₁_dr²c, d²K_dr²c,
     d²W_dr²c, Vc, d²Vdr²c) = taylor_coefficients_low_frequency(
            Kc, Wc, (background, l, ω²))
    H₁₀ = H₁c + 1/2*r₀^2*d²H₁_dr²c
    K₀ = Kc + 1/2*r₀^2*d²K_dr²c
    W₀ = Wc + 1/2*r₀^2*d²W_dr²c
    V₀ = Vc + 1/2*r₀^2*d²Vdr²c

    prob1 = ODEProblem{false}(perturbations_interior_low_frequency,
                              SA[H₁₀, K₀, W₀, V₀], (r₀, rmatch),
                              (background, l, ω²))
    sol1 = solve(prob1, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # Solution 2
    Kc, Wc = zero(ω²), one(ω²)
    (H₁c, d²H₁_dr²c, d²K_dr²c,
     d²W_dr²c, Vc, d²Vdr²c) = taylor_coefficients_low_frequency(
            Kc, Wc, (background, l, ω²))
    H₁₀ = H₁c + 1/2*r₀^2*d²H₁_dr²c
    K₀ = Kc + 1/2*r₀^2*d²K_dr²c
    W₀ = Wc + 1/2*r₀^2*d²W_dr²c
    V₀ = Vc + 1/2*r₀^2*d²Vdr²c

    prob2 = ODEProblem{false}(perturbations_interior_low_frequency,
                              SA[H₁₀, K₀, W₀, V₀], (r₀, rmatch),
                              (background, l, ω²))
    sol2 = solve(prob2, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # from surface
    # Solution 3
    H₁f, Kf, Wf = one(ω²), zero(ω²), zero(ω²)
    Xf = zero(ω²)

    prob3 = ODEProblem{false}(perturbations_interior, SA[H₁f, Kf, Wf, Xf],
                              (R, rmatch), (background, l, ω²))
    sol3 = solve(prob3, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # Solution 4
    H₁f, Kf, Wf = zero(ω²), one(ω²), zero(ω²)
    Xf = zero(ω²)

    prob4 = ODEProblem{false}(perturbations_interior, SA[H₁f, Kf, Wf, Xf],
                              (R, rmatch), (background, l, ω²))
    sol4 = solve(prob4, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # Solution 5
    H₁f, Kf, Wf = zero(ω²), zero(ω²), one(ω²)
    Xf = zero(ω²)

    prob5 = ODEProblem{false}(perturbations_interior, SA[H₁f, Kf, Wf, Xf],
                              (R, rmatch), (background, l, ω²))
    sol5 = solve(prob5, Vern8();
                 abstol=1e-10, reltol=1e-10, save_everystep=save_everystep)

    # solve for coefficients of general solution
    # with W normalised to 1 at surface
    A = SA[0 0 sol3[3, 1] sol4[3, 1] sol5[3, 1]
         sol1[1, end] sol2[1, end] -sol3[1, end] -sol4[1, end] -sol5[1, end]
         sol1[2, end] sol2[2, end] -sol3[2, end] -sol4[2, end] -sol5[2, end]
         sol1[3, end] sol2[3, end] -sol3[3, end] -sol4[3, end] -sol5[3, end]
         X(sol1[end], (background, l, ω²), rmatch) X(sol2[end], (background, l, ω²), rmatch) -sol3[4, end] -sol4[4, end] -sol5[4, end]]
    b = SA[1, 0, 0, 0, 0]
    x = A \ b

    x, sol1, sol2, sol3, sol4, sol5
end

"Return Lagrangian pressure perturbation function."
function X(y, param, r)
    H₁, K, W, V = y
    background, l, ω² = param

    mᵣ = m(background, r)
    pᵣ = p(background, r)
    νᵣ = ν(background, r)
    εᵣ = ε(background.eos, pᵣ)

    expλ = 1 / (1 - 2*mᵣ/r)
    expν = exp(νᵣ)

    dλ_dr = (1 - expλ)/r + 8*π*r*expλ*εᵣ
    dν_dr = (expλ - 1)/r + 8*π*r*expλ*pᵣ

    n = (l + 2)*(l - 1)/2

    H₀ = (
        r^2/expλ*(ω²*r/expν - (n + 1)*dν_dr/2)*H₁
        + (n*r - ω²*r^3/expν - r^2*dν_dr/(4*expλ)*(r*dν_dr - 2))*K
        + 4*π*r^2*(εᵣ + pᵣ)*(
            dν_dr/√expλ*W
            + 2*r*ω²/expν*V
        )
    ) / ((n + 1)*r - r/(2*expλ)*(r*dλ_dr + 2))

    (εᵣ + pᵣ)*(ω²/√expν*V + dν_dr/(2*r)*√(expν/expλ)*W + √expν/2*H₀)
end
