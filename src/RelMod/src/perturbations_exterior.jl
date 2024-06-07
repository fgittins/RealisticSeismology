"""
Functions for mode perturbations in exterior of relativistic star.

Based on method described in Ref. [1]. `transformation` relies on equations in
Refs. [2,3], noting different definitions of metric perturbations.

# Functions
- `V`
- `dV_dr`
- `U`
- `dU_dr`
- `modified_zerilli`
- `integrate_exterior`
- `transformation`
- `Ain`
- `Aout`

# Notes
Assumes geometric units, where G = c = 1. Due to numerical noise at low
frequencies, amplitudes `Ain` and `Aout` are renormalised.

# References
[1] Andersson, Kokkotas and Schutz, "A new numerical approach to the
    oscillation modes of relativistic stars," Mon. Not. R. Astron. Soc. 274
    (4), 1039 (1995).
[2] Detweiler and Lindblom, "On the nonradial pulsations of general
    relativistic stellar models," Astrophys. J. 292, 12 (1985).
[3] Fackerell, "Solutions of Zerilli's equation for even-parity gravitational
    perturbations," Astrophys. J. 166, 197 (1971).
"""

"Effective potential [km^-2] from Zerilli equation."
V(r, M, n) = begin
    2*(1 - 2*M/r)*(n^2*(n + 1)*r^3 + 3*n^2*M*r^2 + 9*n*M^2*r
                   + 9*M^3)/(r^3*(n*r + 3*M)^2)
end

"Derivative of effective potential [km^-3] from Zerilli equation."
dV_dr(r, M, n) = begin 
    (432*M^5 + 54*(10*n - 3)*M^4*r + 18*n*(14*n - 11)*M^3*r^2
     + 6*n^2*(10*n - 13)*M^2*r^3 + 6*n^3*(2*n - 1)*M*r^4
     - 4*n^3*(n + 1)*r^5)/(r^5*(n*r + 3*M)^3)
end

"New effective potential [km^-2]."
U(r, M, n, ω²) = begin
    (1 - 2*M/r)^(-2)*(ω² - V(r, M, n) + 2*M/r^3 - 3*M^2/r^4)
end

"Derivative of new effective potential [km^-3]."
dU_dr(r, M, n, ω²) = begin
    (2*M*(6*M^2 - 8*M*r + 3*r^2 + 2*ω²*r^4) - 4*M*V(r, M, n)*r^4
     + (r - 2*M)*dV_dr(r, M, n)*r^5)/(r*(2*M - r))^3
end

"Exterior perturbation equation."
function modified_zerilli(y, param, ρ)
    q, dq_dρ = y
    R, M, n, ω = param

    θ = -atan(imag(ω)/real(ω))
    r = R + ρ*exp(im*θ)

    d²q_dρ² = 3*dq_dρ^2/(2*q) - 2*q*exp(2*im*θ)*(q^2 - U(r, M, n, ω^2))
    SA[dq_dρ, d²q_dρ²]
end

"Integrate exterior perturbation equation."
function integrate_exterior(R, M, n, ω)
    θ = -atan(imag(ω)/real(ω))
    ρ∞ = 50 / abs(ω)
    r∞ = R + ρ∞*exp(im*θ)
    U∞ = U(r∞, M, n, ω^2)
    dU_dr∞ = dU_dr(r∞, M, n, ω^2)

    prob = ODEProblem{false}(modified_zerilli,
                             SA[√U∞, exp(im*θ)*dU_dr∞/(2*√U∞)], (ρ∞, 0),
                             (R, M, n, ω))
    solve(prob, Vern9(); reltol=1e-15, abstol=1e-15, save_everystep=false)
end

"Transformation from `[H₁, K]` to `[Z, dZdrstar]`."
function transformation(r, x, M, n)
    g = (n*(n + 1)*r^2 + 3*n*M*r + 6*M^2)/(r^2*(n*r + 3*M))
    h = (n*r^2 - 3*n*M*r - 3*M^2)/(r*(r - 2*M)*(n*r + 3*M))
    k = r/(r - 2*M)
    A = SA[1 -k; -g h]/(h - g*k)
    A * x
end

"Amplitude of ingoing radiation [dimensionless]."
Ain(r, q, dq_dr, Z, dZ_drstar, M) = begin
    r*((M/r^2 + (1 - 2*M/r)*(dq_dr/(2*q) + im*q))*Z + dZ_drstar)
end

"Amplitude of outgoing radiation [dimensionless]."
Aout(r, q, dq_dr, Z, dZ_drstar, M) = begin
    -r*((M/r^2 + (1 - 2*M/r)*(dq_dr/(2*q) - im*q))*Z + dZ_drstar)
end
