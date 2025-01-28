"""
Equation-of-state models in β equilibrium at constant entropy per baryon.

# Structs
- `EOS`

# Functions
- `EOS`
- `APR`
- `nb`
- `Base.show`

# Notes
Formulation assumes natural units, where hbar = c = kB = 1.
"""

using Interpolations
using SimpleNonlinearSolve

include("read_table.jl")

mutable struct EOS{IntType₁, IntType₂, IntType₃, IntType₄, IntType₅, IntType₆,
                   IntType₇, IntType₈, IntType₉, IntType₁₀, IntType₁₁,
                   S₁ <: Real, S₂ <: Real, S₃ <: Real}
    const itp₁::IntType₁
    const itp₂::IntType₂
    const itp₃::IntType₃
    const itp₄::IntType₄
    const itp₅::IntType₅
    const itp₆::IntType₆
    const itp₇::IntType₇
    const itp₈::IntType₈
    const itp₉::IntType₉
    const itp₁₀::IntType₁₀
    const itp₁₁::IntType₁₁
    const ymin::S₁
    const ymax::S₂
    const s::S₃     # Entropy per baryon [dimensionless]
    p::Float64      # Pressure [MeV fm^-3]
    nb::Float64     # Baryon-number density [fm^-3]
end

function EOS(s, model, n, m)
    (ss, nbs), (Yₑs, ps, Ts, εs, ∂p_∂Ts, ∂p_∂nbs, ∂p_∂Yₑs,
    ∂s_∂Ts, ∂s_∂nbs, ∂s_∂Yₑs, ∂ε_∂Ts, ∂ε_∂nbs, ∂ε_∂Yₑs) = read_table(
            model, n, m)

    xmin = log(ss[1])
    xmax = log(ss[n])
    δx = (xmax - xmin)/(n - 1)
    xs = xmin:δx:xmax

    ymin = log(nbs[1])
    ymax = log(nbs[m])
    δy = (ymax - ymin)/(m - 1)
    ys = ymin:δy:ymax

    Γ₁s = zeros(n, m)
    ∂logp_∂logTs = zeros(n, m)
    ∂logp_∂lognbs = zeros(n, m)
    ∂logp_∂Yₑs = zeros(n, m)
    ∂logε_∂logTs = zeros(n, m)
    ∂logε_∂lognbs = zeros(n, m)
    ∂logε_∂Yₑs = zeros(n, m)
    for j = 1:m, i = 1:n
        ∂logp_∂logT = Ts[i, j]/ps[i, j]*∂p_∂Ts[i, j]
        ∂logp_∂lognb = nbs[j]/ps[i, j]*∂p_∂nbs[i, j]
        ∂logp_∂Yₑ = 1/ps[i, j]*∂p_∂Yₑs[i, j]

        ∂logs_∂logT = Ts[i, j]/ss[i]*∂s_∂Ts[i, j]
        ∂logs_∂lognb = nbs[j]/ss[i]*∂s_∂nbs[i, j]

        ∂logε_∂logT = Ts[i, j]/εs[i, j]*∂ε_∂Ts[i, j]
        ∂logε_∂lognb = nbs[j]/εs[i, j]*∂ε_∂nbs[i, j]
        ∂logε_∂Yₑ = 1/εs[i, j]*∂ε_∂Yₑs[i, j]

        Γ₁s[i, j] = ∂logp_∂lognb - ∂logp_∂logT*∂logs_∂lognb/∂logs_∂logT
        ∂logp_∂logTs[i, j] = ∂logp_∂logT
        ∂logp_∂lognbs[i, j] = ∂logp_∂lognb
        ∂logp_∂Yₑs[i, j] = ∂logp_∂Yₑ
        ∂logε_∂logTs[i, j] = ∂logε_∂logT
        ∂logε_∂lognbs[i, j] = ∂logε_∂lognb
        ∂logε_∂Yₑs[i, j] = ∂logε_∂Yₑ
    end

    itp₁ = scale(interpolate(log.(ps), BSpline(Cubic())), xs, ys)
    itp₂ = scale(interpolate(log.(εs), BSpline(Cubic())), xs, ys)

    itp₃ = scale(interpolate(Γ₁s, BSpline(Cubic())), xs, ys)

    itp₄ = scale(interpolate(Yₑs, BSpline(Cubic())), xs, ys)
    itp₅ = scale(interpolate(∂logp_∂logTs, BSpline(Linear())), xs, ys)
    itp₆ = scale(interpolate(∂logp_∂lognbs, BSpline(Linear())), xs, ys)
    itp₇ = scale(interpolate(∂logp_∂Yₑs, BSpline(Linear())), xs, ys)
    itp₈ = scale(interpolate(∂logε_∂logTs, BSpline(Linear())), xs, ys)
    itp₉ = scale(interpolate(∂logε_∂lognbs, BSpline(Linear())), xs, ys)
    itp₁₀ = scale(interpolate(∂logε_∂Yₑs, BSpline(Linear())), xs, ys)

    itp₁₁ = scale(interpolate(log.(Ts), BSpline(Linear())), xs, ys)

    p = nb = 0.0

    EOS(itp₁, itp₂, itp₃, itp₄, itp₅, itp₆, itp₇, itp₈, itp₉, itp₁₀, itp₁₁,
        ymin, ymax, s, p, nb)
end

"APR equation of state at constant entropy per baryon."
function APR(s)
    n, m = 21, 220
    model = dirname(Base.active_project()) * "/data/eos/apr_entropy.table"
    EOS(s, model, n, m)    
end

"Baryon-number density [fm^-3] as function of pressure [MeV fm^-3]."
function nb(eos, p)
    prob = IntervalNonlinearProblem{false}(
            (lognb, (eos, logs, logp)) -> eos.itp₁(logs, lognb) - logp,
            (eos.ymin, eos.ymax - 2e-16),
            (eos, log(eos.s), log(p)))
    sol = solve(prob, ITP())
    exp(sol.u)
end

Base.show(io::IO, eos::EOS) = print(io, "EOS with s = $(eos.s)")
