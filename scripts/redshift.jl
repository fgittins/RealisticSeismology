"""
Script calculates f-mode and oscillation spectrum of neutron star with given
baryon mass of Mb = 1.4 Msol, where nuclear matter is described by APR
equation of state at constant redshifted temperature.
"""

using PythonPlot
using Redshift

PythonPlot.rc("text"; usetex=true)
PythonPlot.rc("font"; family="serif", size=16)

include("../src/units.jl")

Mb = 1.4/mass_geometric_to_Msol
eos = APR()
Tᵣ = 0.02

star = calculate_star(eos, Tᵣ, Mb)
println("M = $(star.M*mass_geometric_to_Msol) Msol, R = $(star.R) km")

k₂ = calculate_love_number(star)
println("k₂ = $k₂")

rs = @views star.sol.t[2:end]
νs = @views star.sol[3, 2:end]

Ts = Tᵣ./exp.(νs./2)

fig1, ax1 = subplots()
ax1.plot(log10.(1 .- rs./star.R), log10.(Ts))
ax1.set_xlabel(L"\log_{10}(1 - r / R)")
ax1.set_ylabel(L"$\log_{10}[T(r)$ / MeV$]$")
fig1.tight_layout()
display(fig1)

l = 2
ωguess = (0.075 + 2.9e-5*im)/star.M

ω = solve_eigenfrequency(star, l, ωguess, Simplex())
println("ω M = $(ω*star.M)")
println("Re(ω) / (2π) = $(real(ω)/(2*π)/time_geometric_to_CGS) Hz")

x, sol1, sol2, sol3, sol4, sol5 = Redshift.integrate_interior(
        star, l, ω^2; save_everystep=true)
W(r) = (r ≤ star.R/2
        ? x[1]*sol1(r, idxs=3) + x[2]*sol2(r, idxs=3)
        : x[3]*sol3(r, idxs=3) + x[4]*sol4(r, idxs=3) + x[5]*sol5(r, idxs=3))
rs = LinRange(1e-3, star.R, 200)
ReWs = [real(W(r)) for r ∈ rs]

fig2a, ax2a = subplots()
ax2a.plot(rs, ReWs)
ax2a.set_xlabel(L"$r$ / km")
ax2a.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig2a.tight_layout()
display(fig2a)

fig2b, ax2b = subplots()
ax2b.plot(log10.(1 .- rs./star.R), ReWs)
ax2b.set_xlabel(L"\log_{10}(1 - r / R)")
ax2b.set_ylabel(L"$\mathrm{Re}[W(r)]$ / km$^2$")
fig2b.tight_layout()
display(fig2b)

ωrange = LinRange(0.01, 0.5, 981)./star.M
f(ω) = Redshift.spectrum(star, l, ω)
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
