"""
Script reads in three-parameter APR [1] and DD2 equation-of-state tables [2]
from CompOSE [3-5], interpolates thermodynamic functions and produces
two-parameter tables subject to β equilibrium.

# Structs
- `EOS`

# Functions
- `read_table`
- `APR`
- `DD2`
- `p`
- `s`
- `μΔ`
- `ε`
- `Yₑ`
- `Γ₁`
- `logT`
- `generate_table`
- `generate_table_entropy`
- `make_apr_table`
- `make_dd2_table`
- `make_apr_table_entropy`

# References
[1] Schneider, Constantinou, Muccioli and Prakash,
    "Akmal-Pandharipande-Ravenhall equation of state for simulations of
    supernovae, neutron stars, and binary mergers," Phys. Rev. C 100 (2),
    025803 (2019).
[2] Hempel and Schaffner-Bielich, "A statistical model for a complete supernova
    equation of state," Nucl. Phys. A 837 (3-4), 210 (2010).
[3] Typel, Oertel and Klähn, "CompOSE CompStar online supernova equations of
    state harmonising the concert of nuclear physics and astrophysics
    compose.obspm.fr," Phys. Part. Nucl. 46 (4), 633 (2015).
[4] Oertel, Hempel, Klähn and Typel, "Equations of state for supernovae and
    compact stars," Rev. Mod. Phys. 89 (1), 015007 (2017).
[5] Typel et al., arXiv:2203.03909 [astro-ph.HE].
"""

using DelimitedFiles
using Interpolations
using SimpleNonlinearSolve

"Read three-parameter equation-of-state table from CompOSE."
function read_table(eos)
    Ts = vec(readdlm(eos * "/eos.t"; skipstart=2))::Vector{Float64}
    nbs = vec(readdlm(eos * "/eos.nb"; skipstart=2))::Vector{Float64}
    Yₑs = vec(readdlm(eos * "/eos.yq"; skipstart=2))::Vector{Float64}
    n, m, o = length(Ts), length(nbs), length(Yₑs)

    table, header = readdlm(
            eos * "/eos.thermo"; header=true)::Tuple{Matrix{Float64},
                                                     Matrix{AbstractString}}
    mₙ = parse(Float64, header[1])
    iTs = Int.(@views table[:, 1])
    inbs = Int.(@views table[:, 2])
    iYₑs = Int.(@views table[:, 3])
    ps = zeros(n, m, o)
    ss = zeros(n, m, o)
    μΔs = zeros(n, m, o)
    εs = zeros(n, m, o)
    for i = 1:n*m*o
        ps[iTs[i], inbs[i], iYₑs[i]] = nbs[inbs[i]]*table[i, 4]
        ss[iTs[i], inbs[i], iYₑs[i]] = table[i, 5]
        μΔs[iTs[i], inbs[i], iYₑs[i]] = mₙ*table[i, 8]
        εs[iTs[i], inbs[i], iYₑs[i]] = mₙ*nbs[inbs[i]]*(table[i, 10] + 1)
    end

    (Ts, nbs, Yₑs), (ps, ss, μΔs, εs)
end

"Three-parameter equation of state."
struct EOS{IntType₁, IntType₂, IntType₃, IntType₄,
           S₁ <: Real, S₂ <: Real, S₃ <: Real, S₄ <: Real}
    itp₁::IntType₁
    itp₂::IntType₂
    itp₃::IntType₃
    itp₄::IntType₄
    xmin::S₁
    xmax::S₂
    Yₑmin::S₃
    Yₑmax::S₄
end

"Three-parameter APR equation of state."
function APR()
    (Ts, nbs, Yₑs), (ps, ss, μΔs, εs) = read_table(
            dirname(Base.active_project()) * "/data/eos/apr")

    n, m, o = length(Ts), length(nbs), length(Yₑs)

    xmin = log(Ts[1])
    xmax = log(Ts[n])
    δx = (xmax - xmin)/(n - 1)
    xs = xmin:δx:xmax

    ymin = log(nbs[2])
    ymax = log(nbs[m])
    δy = (ymax - ymin)/(m - 2)
    ys = ymin:δy:ymax

    zmin = Yₑs[1]
    zmax = Yₑs[o]
    δz = (zmax - zmin)/(o - 1)
    zs = zmin:δz:zmax

    itp₁ = scale(interpolate(log.(ps[1:n, 2:m, 1:o]), BSpline(Cubic())),
                 xs, ys, zs)
    itp₂ = scale(interpolate(log.(ss[1:n, 2:m, 1:o]), BSpline(Cubic())),
                 xs, ys, zs)
    itp₃ = scale(interpolate(μΔs[1:n, 2:m, 1:o], BSpline(Cubic())),
                 xs, ys, zs)
    itp₄ = scale(interpolate(log.(εs[1:n, 2:m, 1:o]), BSpline(Cubic())),
                 xs, ys, zs)

    EOS(itp₁, itp₂, itp₃, itp₄, xmin, xmax, zmin, zmax)
end

"Three-parameter DD2 equation of state."
function DD2()
    (Ts, nbs, Yₑs), (ps, ss, μΔs, εs) = read_table(
            dirname(Base.active_project()) * "/data/eos/dd2")

    n, m, o = length(Ts), length(nbs), length(Yₑs)

    xmin = log(Ts[2])
    xmax = log(Ts[n])
    δx = (xmax - xmin)/(n - 2)
    xs = xmin:δx:xmax

    ymin = log(nbs[2])
    ymax = log(nbs[m-2])
    δy = (ymax - ymin)/(m - 4)
    ys = ymin:δy:ymax

    zmin = Yₑs[1]
    zmax = Yₑs[o]
    δz = (zmax - zmin)/(o - 1)
    zs = zmin:δz:zmax

    itp₁ = scale(interpolate(log.(ps[2:n, 2:m-2, 1:o]), BSpline(Cubic())),
                 xs, ys, zs)
    itp₂ = scale(interpolate(log.(ss[2:n, 2:m-2, 1:o]), BSpline(Cubic())),
                 xs, ys, zs)
    itp₃ = scale(interpolate(μΔs[2:n, 2:m-2, 1:o], BSpline(Cubic())),
                 xs, ys, zs)
    itp₄ = scale(interpolate(log.(εs[2:n, 2:m-2, 1:o]), BSpline(Cubic())),
                 xs, ys, zs)

    EOS(itp₁, itp₂, itp₃, itp₄, xmin, xmax, zmin, zmax)
end

"Pressure [MeV fm^-3]."
p(eos, logT, lognb, Yₑ) = exp(eos.itp₁(logT, lognb, Yₑ))
"Entropy per baryon [dimensionless]."
s(eos, logT, lognb, Yₑ) = exp(eos.itp₂(logT, lognb, Yₑ))
"Deviation from β equilibrium [MeV]."
μΔ(eos, logT, lognb, Yₑ) = eos.itp₃(logT, lognb, Yₑ)
"Energy density [MeV fm^-3]."
ε(eos, logT, lognb, Yₑ) = exp(eos.itp₄(logT, lognb, Yₑ))

"Electron fraction [dimensionless] under β equilibrium."
function Yₑ(eos, logT, lognb)
    prob = IntervalNonlinearProblem{false}(
            (Yₑ, (eos, logT, lognb)) -> μΔ(eos, logT, lognb, Yₑ),
            (eos.Yₑmin, eos.Yₑmax), (eos, logT, lognb))
    sol = solve(prob, ITP(); abstol=1e-15, reltol=1e-15)
    sol.u
end

p(eos, logT, lognb) = p(eos, logT, lognb, Yₑ(eos, logT, lognb))
s(eos, logT, lognb) = s(eos, logT, lognb, Yₑ(eos, logT, lognb))
μΔ(eos, logT, lognb) = μΔ(eos, logT, lognb, Yₑ(eos, logT, lognb))
ε(eos, logT, lognb) = ε(eos, logT, lognb, Yₑ(eos, logT, lognb))

"Adiabatic index [dimensionless]."
function Γ₁(eos, logT, lognb, Yₑ)
    ∂logp_∂logT, ∂logp_∂lognb, ∂logp_∂Yₑ = gradient(eos.itp₁, logT, lognb, Yₑ)
    ∂logs_∂logT, ∂logs_∂lognb, ∂logs_∂Yₑ = gradient(eos.itp₂, logT, lognb, Yₑ)
    ∂logp_∂lognb - ∂logp_∂logT*∂logs_∂lognb/∂logs_∂logT
end

Γ₁(eos, logT, lognb) = Γ₁(eos, logT, lognb, Yₑ(eos, logT, lognb))

function logT(eos, logs, lognb)
    prob = IntervalNonlinearProblem{false}(
            (logT, (eos, logs, lognb)) -> log(s(eos, logT, lognb)) - logs,
            (eos.xmin, eos.xmax - 1e-15), (eos, logs, lognb))
    sol = solve(prob, ITP(); abstol=1e-15, reltol=1e-15)
    sol.u
end

"Generate two-parameter equation-of-state table in β equilibrium."
function generate_table(eos, n, m, Tmin, Tmax, nbmin, nbmax)
    xmin = log(Tmin)
    xmax = log(Tmax)
    δx = (xmax - xmin)/(n - 1)
    xs = xmin:δx:xmax

    ymin = log(nbmin)
    ymax = log(nbmax)
    δy = (ymax - ymin)/(m - 1)
    ys = ymin:δy:ymax

    table = zeros(n*m, 15)
    for i = 1:n, j = 1:m
        x, y = xs[i], ys[j]
        z = Yₑ(eos, x, y)   # β equilibrium on each row
        a = p(eos, x, y, z)
        b = s(eos, x, y, z)
        c = ε(eos, x, y, z)

        # T / MeV
        table[m*(i - 1) + j, 1] = exp(x)
        # nb / fm^-3
        table[m*(i - 1) + j, 2] = exp(y)

        # Yₑ
        table[m*(i - 1) + j, 3] = z
        # p / MeV fm^-3
        table[m*(i - 1) + j, 4] = a
        # s
        table[m*(i - 1) + j, 5] = b
        # ε / MeV fm^-3
        table[m*(i - 1) + j, 6] = c

        ∂logp_∂logT, ∂logp_∂lognb, ∂logp_∂Yₑ = gradient(eos.itp₁, x, y, z)
        # ∂p_∂T / fm^-3
        table[m*(i - 1) + j, 7] = a/exp(x)*∂logp_∂logT
        # ∂p_∂nb / MeV
        table[m*(i - 1) + j, 8] = a/exp(y)*∂logp_∂lognb
        # ∂p_∂Yₑ / MeV fm^-3
        table[m*(i - 1) + j, 9] = a*∂logp_∂Yₑ

        ∂logs_∂logT, ∂logs_∂lognb, ∂logs_∂Yₑ = gradient(eos.itp₂, x, y, z)
        # ∂s_∂T / MeV^-1
        table[m*(i - 1) + j, 10] = b/exp(x)*∂logs_∂logT
        # ∂s_∂nb / fm^3
        table[m*(i - 1) + j, 11] = b/exp(y)*∂logs_∂lognb
        # ∂s_∂Yₑ
        table[m*(i - 1) + j, 12] = b*∂logs_∂Yₑ

        ∂logε_∂logT, ∂logε_∂lognb, ∂logε_∂Yₑ = gradient(eos.itp₄, x, y, z)
        # ∂ε_∂T / fm^-3
        table[m*(i - 1) + j, 13] = c/exp(x)*∂logε_∂logT
        # ∂ε_∂nb / MeV
        table[m*(i - 1) + j, 14] = c/exp(y)*∂logε_∂lognb
        # ∂ε_∂Yₑ / MeV fm^-3
        table[m*(i - 1) + j, 15] = c*∂logε_∂Yₑ
    end
    table
end

"""
Generate two-parameter equation-of-state table in β equilibrium with entropy
per baryon as independent variable.
"""
function generate_table_entropy(eos, n, m, smin, smax, nbmin, nbmax)
    Xmin = log(smin)
    Xmax = log(smax)
    δX = (Xmax - Xmin)/(n - 1)
    Xs = Xmin:δX:Xmax

    ymin = log(nbmin)
    ymax = log(nbmax)
    δy = (ymax - ymin)/(m - 1)
    ys = ymin:δy:ymax

    table = zeros(n*m, 15)
    for i = 1:n, j = 1:m
        X, y = Xs[i], ys[j]
        x = logT(eos, X, y)
        z = Yₑ(eos, x, y)   # β equilibrium on each row
        a = p(eos, x, y, z)
        b = exp(X)
        c = ε(eos, x, y, z)

        # s
        table[m*(i - 1) + j, 1] = exp(X)
        # nb / fm^-3
        table[m*(i - 1) + j, 2] = exp(y)

        # Yₑ
        table[m*(i - 1) + j, 3] = z
        # p / MeV fm^-3
        table[m*(i - 1) + j, 4] = a
        # T / MeV
        table[m*(i - 1) + j, 5] = exp(x)
        # ε / MeV fm^-3
        table[m*(i - 1) + j, 6] = c

        ∂logp_∂logT, ∂logp_∂lognb, ∂logp_∂Yₑ = gradient(eos.itp₁, x, y, z)
        # ∂p_∂T / fm^-3
        table[m*(i - 1) + j, 7] = a/exp(x)*∂logp_∂logT
        # ∂p_∂nb / MeV
        table[m*(i - 1) + j, 8] = a/exp(y)*∂logp_∂lognb
        # ∂p_∂Yₑ / MeV fm^-3
        table[m*(i - 1) + j, 9] = a*∂logp_∂Yₑ

        ∂logs_∂logT, ∂logs_∂lognb, ∂logs_∂Yₑ = gradient(eos.itp₂, x, y, z)
        # ∂s_∂T / MeV^-1
        table[m*(i - 1) + j, 10] = b/exp(x)*∂logs_∂logT
        # ∂s_∂nb / fm^3
        table[m*(i - 1) + j, 11] = b/exp(y)*∂logs_∂lognb
        # ∂s_∂Yₑ
        table[m*(i - 1) + j, 12] = b*∂logs_∂Yₑ

        ∂logε_∂logT, ∂logε_∂lognb, ∂logε_∂Yₑ = gradient(eos.itp₄, x, y, z)
        # ∂ε_∂T / fm^-3
        table[m*(i - 1) + j, 13] = c/exp(x)*∂logε_∂logT
        # ∂ε_∂nb / MeV
        table[m*(i - 1) + j, 14] = c/exp(y)*∂logε_∂lognb
        # ∂ε_∂Yₑ / MeV fm^-3
        table[m*(i - 1) + j, 15] = c*∂logε_∂Yₑ
    end
    table
end

function make_apr_table()
    apr = APR()

    # range to enforce β equilibrium
    n, m = 132, 220
    Tmin = 9.9770006382255347e-3
    Tmax = 232.09547139699004
    nbmin = 1.2618275141636872e-7
    nbmax = 2.5176768870016693

    table = generate_table(apr, n, m, Tmin, Tmax, nbmin, nbmax)

    open(dirname(Base.active_project()) * "/data/eos/apr.table", "w") do io
        writedlm(io, table)
    end

    nothing
end

function make_dd2_table()
    dd2 = DD2()

    # range to enforce β equilibrium
    n, m = 80, 323
    Tmin = 1.0964782e-1
    Tmax = 1.5848932e2
    nbmin = 1.0964782e-12
    nbmax = 8.3176377

    table = generate_table(dd2, n, m, Tmin, Tmax, nbmin, nbmax)

    open(dirname(Base.active_project()) * "/data/eos/dd2.table", "w") do io
        writedlm(io, table)
    end

    nothing
end

function make_apr_table_entropy()
    apr = APR()

    # range to enforce β equilibrium
    n, m = 21, 220
    smin = 0.1
    smax = 5.0
    nbmin = 1.2618275141636872e-7
    nbmax = 2.5176768870016693

    table = generate_table_entropy(apr, n, m, smin, smax, nbmin, nbmax)

    open(dirname(Base.active_project()) * "/data/eos/apr_entropy.table", "w") do io
        writedlm(io, table)
    end

    nothing
end
