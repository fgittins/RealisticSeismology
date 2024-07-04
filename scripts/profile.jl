"""
Script calculates f-mode and oscillation spectrum of neutron star with given
baryon mass of Mb = 1.4 Msol, where nuclear matter is described by APR equation
of state with temperature profile from numerical-relativity simulation.
"""

using PythonPlot
using RelMod

PythonPlot.rc("text"; usetex=true)
PythonPlot.rc("font"; family="serif", size=16)

include("../src/baryon_mass.jl")
include("../src/love.jl")
include("../src/sim.jl")
include("../src/units.jl")

function RelMod.ε(sim::Sim, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ sim.p
        sim.p = pnatural
        sim.nb = nb(sim, pnatural)
        sim.T = sim.itp₁₁(log(sim.nb))
    end
    εnatural = exp(sim.itp₂(log(sim.T), log(sim.nb)))
    εnatural/pressure_geometric_to_natural
end

function RelMod.Γ(sim::Sim, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ sim.p
        sim.p = pnatural
        sim.nb = nb(sim, pnatural)
        sim.T = sim.itp₁₁(log(sim.nb))
    end
    εnatural = exp(sim.itp₂(log(sim.T), log(sim.nb)))

    ∂Yₑ_∂logT, ∂Yₑ_∂lognb = gradient(sim.itp₄, log(sim.T), log(sim.nb))

    dT_dlognb = gradient(sim.itp₁₁, log(sim.nb))[1]
    dlogT_dlognb = 1/sim.T*dT_dlognb

    ∂logp_∂logT = sim.itp₅(log(sim.T), log(sim.nb))
    ∂logp_∂lognb = sim.itp₆(log(sim.T), log(sim.nb))
    ∂logp_∂Yₑ = sim.itp₇(log(sim.T), log(sim.nb))
    dlogp_dlognb = (∂logp_∂lognb + ∂logp_∂Yₑ*∂Yₑ_∂lognb
                    + (∂logp_∂logT + ∂logp_∂Yₑ*∂Yₑ_∂logT)*dlogT_dlognb)

    ∂logε_∂logT = sim.itp₈(log(sim.T), log(sim.nb))
    ∂logε_∂lognb = sim.itp₉(log(sim.T), log(sim.nb))
    ∂logε_∂Yₑ = sim.itp₁₀(log(sim.T), log(sim.nb))
    dlogε_dlognb = (∂logε_∂lognb + ∂logε_∂Yₑ*∂Yₑ_∂lognb
                    + (∂logε_∂logT + ∂logε_∂Yₑ*∂Yₑ_∂logT)*dlogT_dlognb)

    (εnatural + sim.p)/εnatural*dlogp_dlognb/dlogε_dlognb
end

function RelMod.Γ₁(sim::Sim, p)
    pnatural = pressure_geometric_to_natural*p
    if pnatural ≠ sim.p
        sim.p = pnatural
        sim.nb = nb(sim, pnatural)
        sim.T = sim.itp₁₁(log(sim.nb))
    end
    sim.itp₃(log(sim.T), log(sim.nb))
end

Mb = 1.4/mass_geometric_to_Msol
eos = APR()
pf = exp(eos.itp₁(log(eos.itp₁₁(eos.ymin)), eos.ymin))/pressure_geometric_to_natural

star = calculate_star(eos, Mb, pf)
println("M = $(star.M*mass_geometric_to_Msol) Msol, R = $(star.R) km")

k₂ = calculate_love_number(star)
println("k₂ = $k₂")

l = 2
ωguess = (0.063 + 1.2e-7*im)/star.M

ω = solve_eigenfrequency(star, l, ωguess, Simplex())
println("ω M = $(ω*star.M)")
println("Re(ω) / (2π) = $(real(ω)/(2*π)/time_geometric_to_CGS) Hz")

x, sol1, sol2, sol3, sol4, sol5 = RelMod.integrate_interior(
        star, l, ω^2; save_everystep=true)
W(r) = (r ≤ star.R/2
        ? x[1]*sol1(r, idxs=3) + x[2]*sol2(r, idxs=3)
        : x[3]*sol3(r, idxs=3) + x[4]*sol4(r, idxs=3) + x[5]*sol5(r, idxs=3))
rs = LinRange(1e-3, star.R, 200)
ReWs = [real(W(r)) for r ∈ rs]

fig1, ax1 = subplots()
ax1.plot(rs, ReWs)
ax1.set_xlabel(L"$r$ / km")
ax1.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig1.tight_layout()
display(fig1)

ωrange = LinRange(0.01, 0.5, 981)./star.M
f(ω) = RelMod.spectrum(star, l, ω)
frange = f.(ωrange)

fig2a, ax2a = subplots()
ax2a.plot(ωrange.*star.M, log10.(abs.(frange)))
ax2a.set_xlabel(L"\omega M")
ax2a.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig2a.tight_layout()
display(fig2a)

fig2b, ax2b = subplots()
ax2b.plot(ωrange./(2*π*time_geometric_to_CGS)./1e3, log10.(abs.(frange)))
ax2b.set_xlabel(L"$\omega / (2 \pi)$ / kHz")
ax2b.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig2b.tight_layout()
display(fig2b)

ωrange = LinRange(0.001, 0.025, 2401)./star.M
frange = f.(ωrange)

fig3a, ax3a = subplots()
ax3a.plot(ωrange.*star.M, log10.(abs.(frange)))
ax3a.set_xlabel(L"\omega M")
ax3a.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig3a.tight_layout()
display(fig3a)

fig3b, ax3b = subplots()
ax3b.plot(ωrange./(2*π*time_geometric_to_CGS), log10.(abs.(frange)))
ax3b.set_xlabel(L"$\omega / (2 \pi)$ / Hz")
ax3b.set_ylabel(L"\log_{10}(|\tilde{A}_\mathrm{in}|)")
fig3b.tight_layout()
display(fig3b)

ωrange = LinRange(0.0008, 0.008, 3601)./star.M
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
