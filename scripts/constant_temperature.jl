"""
Script calculates f-mode and oscillation spectrum of neutron star with given
baryon mass of Mb = 1.4 Msol, where nuclear matter is described by APR
equation of state at constant temperature.
"""

using PythonPlot
using RelMod

PythonPlot.rc("text"; usetex=true)
PythonPlot.rc("font"; family="serif", size=16)

include("../src/eos.jl")
include("../src/baryon_mass.jl")
include("../src/love.jl")
include("../src/units.jl")

function RelMod.ε(eos::EOS, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ eos.p
        eos.p = pnatural
        eos.nb = nb(eos, pnatural)
    end
    εnatural = exp(eos.itp₂(log(eos.T), log(eos.nb)))
    εnatural/pressure_geometric_to_natural
end

function RelMod.Γ(eos::EOS, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ eos.p
        eos.p = pnatural
        eos.nb = nb(eos, pnatural)
    end
    εnatural = exp(eos.itp₂(log(eos.T), log(eos.nb)))

    ∂Yₑ_∂logT, ∂Yₑ_∂lognb = gradient(eos.itp₄, log(eos.T), log(eos.nb))

    ∂logp_∂lognb = eos.itp₆(log(eos.T), log(eos.nb))
    ∂logp_∂Yₑ = eos.itp₇(log(eos.T), log(eos.nb))
    dlogp_dlognb = ∂logp_∂lognb + ∂logp_∂Yₑ*∂Yₑ_∂lognb

    ∂logε_∂lognb = eos.itp₉(log(eos.T), log(eos.nb))
    ∂logε_∂Yₑ = eos.itp₁₀(log(eos.T), log(eos.nb))
    dlogε_dlognb = ∂logε_∂lognb + ∂logε_∂Yₑ*∂Yₑ_∂lognb

    (εnatural + eos.p)/εnatural*dlogp_dlognb/dlogε_dlognb
end

function RelMod.Γ₁(eos::EOS, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ eos.p
        eos.p = pnatural
        eos.nb = nb(eos, pnatural)
    end
    eos.itp₃(log(eos.T), log(eos.nb))
end

T = 0.02
Mb = 1.4/mass_geometric_to_Msol
eos = APR(T)
pf = exp(eos.itp₁(log(eos.T), eos.ymin))/pressure_geometric_to_natural

star = calculate_star(eos, Mb, pf)
println("M = $(star.M*mass_geometric_to_Msol) Msol, R = $(star.R) km")

k₂ = calculate_love_number(star)
println("k₂ = $k₂")

rs = @views star.sol.t[2:end]
ps = @views star.sol[2, 2:end]

ε(p) = RelMod.ε(eos, p)
εs = ε.(ps)

fig1, ax1 = subplots()
ax1.plot(log10.(1 .- rs./star.R), log10.(εs.*pressure_geometric_to_natural))
ax1.axvline(log10(1 - rs[26]/star.R);
            alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax1.axvline(log10(1 - rs[101]/star.R);
            alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax1.set_xlabel(L"\log_{10}(1 - r / R)")
ax1.set_ylabel(L"$\log_{10}[\varepsilon(r)$ / MeV fm$^{-3}]$")
fig1.tight_layout()
display(fig1)

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

fig2a, ax2a = subplots()
ax2a.plot(rs, ReWs)
ax2a.set_xlabel(L"$r$ / km")
ax2a.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig2a.tight_layout()
display(fig2a)

fig2b, ax2b = subplots()
ax2b.plot(log10.(1 .- rs./star.R), ReWs)
ax2b.axvline(log10(1 - rs[26]/star.R);
             alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax2b.axvline(log10(1 - rs[101]/star.R);
             alpha=0.5, color="tab:gray", linestyle="--", linewidth=1)
ax2b.set_xlabel(L"\log_{10}(1 - r / R)")
ax2b.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig2b.tight_layout()
display(fig2b)

ωrange = LinRange(0.01, 0.5, 981)./star.M
f(ω) = RelMod.spectrum(star, l, ω)
frange = f.(ωrange)

fig3a, ax3a = subplots()
ax3a.plot(ωrange.*star.M, log10.(abs.(frange)))
ax3a.set_xlabel(L"\omega M")
ax3a.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig3a.tight_layout()
display(fig3a)

fig3b, ax3b = subplots()
ax3b.plot(ωrange./(2*π*time_geometric_to_CGS)./1e3, log10.(abs.(frange)))
ax3b.set_xlabel(L"$\omega / (2 \pi)$ / kHz")
ax3b.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig3b.tight_layout()
display(fig3b)

ωrange = LinRange(0.001, 0.025, 2401)./star.M
frange = f.(ωrange)

fig4a, ax4a = subplots()
ax4a.plot(ωrange.*star.M, log10.(abs.(frange)))
ax4a.set_xlabel(L"\omega M")
ax4a.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig4a.tight_layout()
display(fig4a)

fig4b, ax4b = subplots()
ax4b.plot(ωrange./(2*π*time_geometric_to_CGS), log10.(abs.(frange)))
ax4b.set_xlabel(L"$\omega / (2 \pi)$ / Hz")
ax4b.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig4b.tight_layout()
display(fig4b)

ωrange = LinRange(0.0008, 0.008, 3601)./star.M
frange = f.(ωrange)

fig5a, ax5a = subplots()
ax5a.plot(ωrange.*star.M, log10.(abs.(frange)))
ax5a.set_xlabel(L"\omega M")
ax5a.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig5a.tight_layout()
display(fig5a)

fig5b, ax5b = subplots()
ax5b.plot(ωrange./(2*π*time_geometric_to_CGS), log10.(abs.(frange)))
ax5b.set_xlabel(L"$\omega / (2 \pi)$ / Hz")
ax5b.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig5b.tight_layout()
display(fig5b)
