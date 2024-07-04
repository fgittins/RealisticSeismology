"""
Equation of state with temperature profile from simulation.

# Structs
- `Sim`

# Functions
- `Sim`
- `APR`
- `DD2`
- `nb`

# Notes
Model is taken from CompOSE repository, which provides data in tabulated form.
Table is initially three-dimensional in temperature `T` [MeV], baryon-number
density `nb` [fm^-3] and electron fraction `Yₑ` [dimensionless].

Table is put in two-dimensional form by enforcing β equilibrium, solving for
`Yₑ` at each `(T, nb)`. Then, model is interpolated in these two dimensions.

This implementation takes profile of `T = T[log10(nb/nbsat)]` from
numerical-relativity simulation. Baryon-number density is determined by
inverting `p = p(nb)`.

Formulation assumes natural units, where hbar = c = kB = 1.
"""

using Interpolations
using SimpleNonlinearSolve

include("read_table.jl")

"Equation of state with temperature profile from simulation."
mutable struct Sim{IntType₁, IntType₂, IntType₃, IntType₄, IntType₅, IntType₆,
                   IntType₇, IntType₈, IntType₉, IntType₁₀,
                   T₁ <: Real, T₂ <: Real,
                   IntType₁₁}
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
    const ymin::T₁
    const ymax::T₂
    T::Float64              # Temperature [MeV]
    p::Float64              # Pressure [MeV fm^-3]
    nb::Float64             # Baryon-number density [fm^-3]
    const itp₁₁::IntType₁₁  # Temperature profile [MeV]
end

function Sim(n, m, model, sim)
    (Ts, nbs), (Yₑs, ps, ss, εs, ∂p_∂Ts, ∂p_∂nbs, ∂p_∂Yₑs,
    ∂s_∂Ts, ∂s_∂nbs, ∂s_∂Yₑs, ∂ε_∂Ts, ∂ε_∂nbs, ∂ε_∂Yₑs) = read_table(
            model, n, m)

    xmin = log(Ts[1])
    xmax = log(Ts[n])
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
        ∂logp_∂logT = Ts[i]/ps[i, j]*∂p_∂Ts[i, j]
        ∂logp_∂lognb = nbs[j]/ps[i, j]*∂p_∂nbs[i, j]
        ∂logp_∂Yₑ = 1/ps[i, j]*∂p_∂Yₑs[i, j]

        ∂logs_∂logT = Ts[i]/ss[i, j]*∂s_∂Ts[i, j]
        ∂logs_∂lognb = nbs[j]/ss[i, j]*∂s_∂nbs[i, j]

        ∂logε_∂logT = Ts[i]/εs[i, j]*∂ε_∂Ts[i, j]
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

    T = p = nb = 0.0

    profile = readdlm(sim; comments=true)::Matrix{Float64}

    xs = @view profile[:, 1]
    @. xs = xs*log(10)          # convert to natural logarithm
    nbsat = 0.160               # saturation density [fm^-3]
    lognbs = xs .+ log(nbsat)
    Ts = @view profile[:, 2]
    itp₁₁ = interpolate((lognbs,), Ts, Gridded(Linear()))

    ymin = lognbs[1] ≥ ymin ? lognbs[1] : ymin
    ymax = lognbs[end] ≤ ymax ? lognbs[end] : ymax

    Sim(itp₁, itp₂, itp₃, itp₄, itp₅, itp₆, itp₇, itp₈, itp₉, itp₁₀,
        ymin, ymax, T, p, nb, itp₁₁)
end

"APR equation of state with temperature profile from simulation."
function APR()
    n, m = 132, 220
    model = dirname(Base.active_project()) * "/data/eos/apr.table"
    sim = dirname(Base.active_project()) * "/data/sim/apr.dat"
    Sim(n, m, model, sim)
end

"DD2 equation of state with temperature profile from simulation."
function DD2()
    n, m = 80, 323
    model = dirname(Base.active_project()) * "/data/eos/dd2.table"
    sim = dirname(Base.active_project()) * "/data/sim/dd2.dat"
    Sim(n, m, model, sim)
end

"Baryon-number density [fm^-3] as function of pressure [MeV fm^-3]."
function nb(sim, p)
    prob = IntervalNonlinearProblem{false}(
            (lognb, (sim, logp)) -> sim.itp₁(log(sim.itp₁₁(lognb)), lognb) - logp,
            (sim.ymin, sim.ymax - 2e-16),
            (sim, log(p)))
    sol = solve(prob, ITP())
    exp(sol.u)
end
