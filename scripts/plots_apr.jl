"""
Script reads in two-parameter APR equation-of-state table in β equilibrium and
plots thermodynamic functions.
"""

using Interpolations
using PythonPlot

PythonPlot.rc("text"; usetex=true)
PythonPlot.rc("font"; family="serif", size=16)
PythonPlot.rc("figure"; dpi=140)

include("../src/read_table.jl")

n, m = 132, 220

(Ts, nbs), (Yₑs, ps, ss, εs, ∂p_∂Ts, ∂p_∂nbs, ∂p_∂Yₑs,
∂s_∂Ts, ∂s_∂nbs, ∂s_∂Yₑs, ∂ε_∂Ts, ∂ε_∂nbs, ∂ε_∂Yₑs) = read_table(
        dirname(Base.active_project()) * "/data/eos/apr.table", n, m)

xmin = log(Ts[1])
xmax = log(Ts[end])
δx = (xmax - xmin)/(n - 1)
xs = xmin:δx:xmax

ymin = log(nbs[1])
ymax = log(nbs[end])
δy = (ymax - ymin)/(m - 1)
ys = ymin:δy:ymax

itp = scale(interpolate(Yₑs, BSpline(Cubic())), xs, ys)

Γ₁s = zeros(n, m)
Γs = zeros(n, m)
for j = 1:m, i = 1:n
    ∂Yₑ_∂logT, ∂Yₑ_∂lognb = gradient(itp, xs[i], ys[j])

    ∂logp_∂logT = Ts[i]/ps[i, j]*∂p_∂Ts[i, j]
    ∂logp_∂lognb = nbs[j]/ps[i, j]*∂p_∂nbs[i, j]
    ∂logp_∂Yₑ = 1/ps[i, j]*∂p_∂Yₑs[i, j]
    dlogp_dlognb = ∂logp_∂lognb + ∂logp_∂Yₑ*∂Yₑ_∂lognb

    ∂logs_∂logT = Ts[i]/ss[i, j]*∂s_∂Ts[i, j]
    ∂logs_∂lognb = nbs[j]/ss[i, j]*∂s_∂nbs[i, j]

    ∂logε_∂lognb = nbs[j]/εs[i, j]*∂ε_∂nbs[i, j]
    ∂logε_∂Yₑ = 1/εs[i, j]*∂ε_∂Yₑs[i, j]
    dlogε_dlognb = ∂logε_∂lognb + ∂logε_∂Yₑ*∂Yₑ_∂lognb

    Γ₁s[i, j] = ∂logp_∂lognb - ∂logp_∂logT*∂logs_∂lognb/∂logs_∂logT
    Γs[i, j] = (εs[i, j] + ps[i, j])/εs[i, j]*dlogp_dlognb/dlogε_dlognb
end

fig1, ax1 = subplots()
cs1 = ax1.contourf(log10.(Ts), log10.(nbs), Yₑs')
ax1.set_xlabel(L"$\log_{10}(T$ / MeV$)$")
ax1.set_ylabel(L"$\log_{10}(n_\mathrm{b}$ / fm$^{-3})$")
cbar1 = fig1.colorbar(cs1)
cbar1.ax.set_ylabel(L"Y_\mathrm{e}")
fig1.tight_layout()
display(fig1)

fig2, ax2 = subplots()
cs2 = ax2.contourf(log10.(Ts), log10.(nbs), log10.(ps'))
ax2.set_xlabel(L"$\log_{10}(T$ / MeV$)$")
ax2.set_ylabel(L"$\log_{10}(n_\mathrm{b}$ / fm$^{-3})$")
cbar2 = fig2.colorbar(cs2)
cbar2.ax.set_ylabel(L"$\log_{10}(p$ / MeV fm$^{-3})$")
fig2.tight_layout()
display(fig2)

fig3, ax3 = subplots()
cs3 = ax3.contourf(log10.(Ts), log10.(nbs), log10.(εs'))
ax3.set_xlabel(L"$\log_{10}(T$ / MeV$)$")
ax3.set_ylabel(L"$\log_{10}(n_\mathrm{b}$ / fm$^{-3})$")
cbar3 = fig3.colorbar(cs3)
cbar3.ax.set_ylabel(L"$\log_{10}(\varepsilon$ / MeV fm$^{-3})$")
fig3.tight_layout()
display(fig3)

fig4, ax4 = subplots()
cs4 = ax4.contourf(log10.(Ts), log10.(nbs), Γs')
ax4.set_xlabel(L"$\log_{10}(T$ / MeV$)$")
ax4.set_ylabel(L"$\log_{10}(n_\mathrm{b}$ / fm$^{-3})$")
cbar4 = fig4.colorbar(cs4)
cbar4.ax.set_ylabel(L"\Gamma")
ax4.set_title(L"T = \mathrm{const}")
fig4.tight_layout()
display(fig4)

fig5, ax5 = subplots()
cs5 = ax5.contourf(log10.(Ts), log10.(nbs), Γ₁s')
ax5.set_xlabel(L"$\log_{10}(T$ / MeV$)$")
ax5.set_ylabel(L"$\log_{10}(n_\mathrm{b}$ / fm$^{-3})$")
cbar5 = fig5.colorbar(cs5)
cbar5.ax.set_ylabel(L"\Gamma_1")
fig5.tight_layout()
display(fig5)

fig6, ax6 = subplots()
ax6.plot(log10.(nbs), Γs[1, :]; label=L"\Gamma")
ax6.plot(log10.(nbs), Γ₁s[1, :]; label=L"\Gamma_1")
ax6.set_xlabel(L"$\log_{10}(n_\mathrm{b}$ / fm$^{-3})$")
ax6.legend()
ax6.set_title("\$T = $(Ts[1]) \\, \\mathrm{MeV} = \\mathrm{const}\$")
fig6.tight_layout()
display(fig6)
