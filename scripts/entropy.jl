"""
Script calculates f-mode neutron star with given baryon mass of Mb = 1.4 Msol, 
where nuclear matter is described by APR equation of state at constant entropy
per baryon.
"""

using PythonPlot
using RelMod

PythonPlot.rc("text"; usetex=true)
PythonPlot.rc("font"; family="serif", size=16)

include("../src/eos_entropy.jl")
include("../src/baryon_mass.jl")
include("../src/love.jl")
include("../src/units.jl")

function RelMod.ε(eos::EOS, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ eos.p
        eos.p = pnatural
        eos.nb = nb(eos, pnatural)
    end
    εnatural = exp(eos.itp₂(log(eos.s), log(eos.nb)))
    εnatural/pressure_geometric_to_natural
end

function RelMod.Γ(eos::EOS, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ eos.p
        eos.p = pnatural
        eos.nb = nb(eos, pnatural)
    end
    εnatural = exp(eos.itp₂(log(eos.s), log(eos.nb)))

    ∂Yₑ_∂logT, ∂Yₑ_∂lognb = gradient(eos.itp₄, log(eos.s), log(eos.nb))

    ∂logp_∂lognb = eos.itp₆(log(eos.s), log(eos.nb))
    ∂logp_∂Yₑ = eos.itp₇(log(eos.s), log(eos.nb))
    dlogp_dlognb = ∂logp_∂lognb + ∂logp_∂Yₑ*∂Yₑ_∂lognb

    ∂logε_∂lognb = eos.itp₉(log(eos.s), log(eos.nb))
    ∂logε_∂Yₑ = eos.itp₁₀(log(eos.s), log(eos.nb))
    dlogε_dlognb = ∂logε_∂lognb + ∂logε_∂Yₑ*∂Yₑ_∂lognb

    (εnatural + eos.p)/εnatural*dlogp_dlognb/dlogε_dlognb
end

function RelMod.Γ₁(eos::EOS, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ eos.p
        eos.p = pnatural
        eos.nb = nb(eos, pnatural)
    end
    eos.itp₃(log(eos.s), log(eos.nb))
end

T(eos::EOS, nb) = exp(eos.itp₁₁(log(eos.s), log(nb)))

s = 0.3
Mb = 1.4/mass_geometric_to_Msol
eos = APR(s)
pf = exp(eos.itp₁(log(eos.s), eos.ymin))/pressure_geometric_to_natural

star = calculate_star(eos, Mb, pf)
println("M = $(star.M*mass_geometric_to_Msol) Msol, R = $(star.R) km")

k₂ = calculate_love_number(star)
println("k₂ = $k₂")

rs = @views star.sol.t[2:end]
ps = @views star.sol[2, 2:end]

nb(p) = nb(eos, pressure_geometric_to_natural*p)
T(nb) = T(eos, nb)
nbs = nb.(ps)
Ts = T.(nbs)

fig1, ax1 = subplots()
ax1.plot(log10.(1 .- rs./star.R), log10.(nbs))
ax1.axvline(log10(1 - rs[34]/star.R);
            alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax1.axvline(log10(1 - rs[92]/star.R);
            alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax1.set_xlabel(L"\log_{10}(1 - r / R)")
ax1.set_ylabel(L"$\log_{10}[n_\mathrm{b}(r)$ / fm$^{-3}]$")
fig1.tight_layout()
display(fig1)

fig2, ax2 = subplots()
ax2.plot(log10.(1 .- rs./star.R), log10.(Ts))
ax2.axvline(log10(1 - rs[34]/star.R);
            alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax2.axvline(log10(1 - rs[92]/star.R);
            alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax2.set_xlabel(L"\log_{10}(1 - r / R)")
ax2.set_ylabel(L"$\log_{10}[T(r)$ / MeV$]$")
fig2.tight_layout()
display(fig2)

l = 2
ωguess = (0.075 + 2.9e-5*im)/star.M

ω = solve_eigenfrequency(star, l, ωguess, Simplex())
println("ω M = $(ω*star.M)")
println("Re(ω) / (2π) = $(real(ω)/(2*π)/time_geometric_to_CGS) Hz")

x, sol1, sol2, sol3, sol4, sol5 = RelMod.integrate_interior(
        star, l, ω^2; save_everystep=true)
W(r) = (r ≤ star.R/2
        ? x[1]*sol1(r, idxs=3) + x[2]*sol2(r, idxs=3)
        : x[3]*sol3(r, idxs=3) + x[4]*sol4(r, idxs=3) + x[5]*sol5(r, idxs=3))
ReWs = [real(W(r)) for r ∈ rs]

fig3a, ax3a = subplots()
ax3a.plot(rs, ReWs)
ax3a.set_xlabel(L"$r$ / km")
ax3a.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig3a.tight_layout()
display(fig3a)

fig3b, ax3b = subplots()
ax3b.plot(log10.(1 .- rs./star.R), ReWs)
ax3b.axvline(log10(1 - rs[34]/star.R);
             alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax3b.axvline(log10(1 - rs[92]/star.R);
             alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax3b.set_xlabel(L"\log_{10}(1 - r / R)")
ax3b.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig3b.tight_layout()
display(fig3b)
