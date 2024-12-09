"""
Two-parameter APR equation of state.

# Structs
- `APR`

# Functions
- `APR`
- `nb`
- `ε`
- `Γ`
- `Γ₁`
- `Base.show`
"""

mutable struct APR{IntType₁, IntType₂, IntType₃, IntType₄, IntType₅, IntType₆,
                   IntType₇, IntType₈, IntType₉, IntType₁₀,
                   S₁ <: Real, S₂ <: Real}
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
    const ymin::S₁
    const ymax::S₂
    T::Float64      # Temperature [MeV]
    p::Float64      # Pressure [MeV fm^-3]
    nb::Float64     # Baryon-number density [fm^-3]
end

"Two-parameter APR equation of state."
function APR()
    n, m = 132, 220

    (Ts, nbs), (Yₑs, ps, ss, εs, ∂p_∂Ts, ∂p_∂nbs, ∂p_∂Yₑs,
    ∂s_∂Ts, ∂s_∂nbs, ∂s_∂Yₑs, ∂ε_∂Ts, ∂ε_∂nbs, ∂ε_∂Yₑs) = read_table(
            dirname(Base.active_project()) * "/data/eos/apr.table", n, m)

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

    APR(itp₁, itp₂, itp₃, itp₄, itp₅, itp₆, itp₇, itp₈, itp₉, itp₁₀,
        ymin, ymax, T, p, nb)
end

"""
Baryon-number density [fm^-3] as function of temperature [MeV] and pressure
[MeV fm^-3].
"""
function nb(apr, T, p)
    prob = IntervalNonlinearProblem{false}(
            (lognb, (apr, logT, logp)) -> apr.itp₁(logT, lognb) - logp,
            (apr.ymin, apr.ymax - 2e-16),
            (apr, log(T), log(p)))
    sol = solve(prob, ITP())
    exp(sol.u)
end

"Energy density [km^-2] as function of temperature [MeV] and pressure [km^-2]."
function ε(apr, T, p)
    pnatural = pressure_geometric_to_natural*p
    if T ≠ apr.T && pnatural ≠ apr.p
        apr.T = T
        apr.p = pnatural
        apr.nb = nb(apr, T, pnatural)
    end
    εnatural = exp(apr.itp₂(log(apr.T), log(apr.nb)))
    εnatural/pressure_geometric_to_natural
end

"""
Background index [dimensionless] for constant redshifted temperature as
function of temperature [MeV] and pressure [km^-2].
"""
function Γ(apr, T, p)
    pnatural = pressure_geometric_to_natural*p
    if T ≠ apr.T && pnatural ≠ apr.p
        apr.T = T
        apr.p = pnatural
        apr.nb = nb(apr, T, pnatural)
    end
    εnatural = exp(apr.itp₂(log(apr.T), log(apr.nb)))

    ∂Yₑ_∂logT, ∂Yₑ_∂lognb = gradient(apr.itp₄, log(apr.T), log(apr.nb))

    ∂logp_∂logT = apr.itp₅(log(apr.T), log(apr.nb))
    ∂logp_∂lognb = apr.itp₆(log(apr.T), log(apr.nb))
    ∂logp_∂Yₑ = apr.itp₇(log(apr.T), log(apr.nb))
    dlogp_dlognb = (∂logp_∂lognb + ∂logp_∂Yₑ*∂Yₑ_∂lognb)/(1 
            - apr.p/(εnatural + apr.p)*(∂logp_∂logT + ∂logp_∂Yₑ*∂Yₑ_∂logT))

    dlogT_dlognb = apr.p/(εnatural + apr.p)*dlogp_dlognb

    ∂logε_∂logT = apr.itp₈(log(apr.T), log(apr.nb))
    ∂logε_∂lognb = apr.itp₉(log(apr.T), log(apr.nb))
    ∂logε_∂Yₑ = apr.itp₁₀(log(apr.T), log(apr.nb))
    dlogε_dlognb = (∂logε_∂lognb + ∂logε_∂Yₑ*∂Yₑ_∂lognb
                    + (∂logε_∂logT + ∂logε_∂Yₑ*∂Yₑ_∂logT)*dlogT_dlognb)

    (εnatural + apr.p)/εnatural*dlogp_dlognb/dlogε_dlognb
end

"""
Adiabatic index [dimensionless] as function of temperature [MeV] and pressure
[km^-2].
"""
function Γ₁(apr, T, p)
    pnatural = pressure_geometric_to_natural*p
    if T ≠ apr.T && pnatural ≠ apr.p
        apr.T = T
        apr.p = pnatural
        apr.nb = nb(apr, T, pnatural)
    end
    apr.itp₃(log(apr.T), log(apr.nb))
end

Base.show(io::IO, apr::APR) = print(io, "APR")
