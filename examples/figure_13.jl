using SpecialFunctions: erfcx
import RationalApproximations: locate_extrema_idx, RationalFunc, aaa_v2, minimax
using PyPlot

G(t) = (1-t)/(1+t)
α = 0.5
β = 1.0
fα(t) = erfcx(G(t))

t = range(-1, 1, length=800)
max_m = 5
tol = 1e-10
F = fα.(t)
r, err = aaa_v2(F, t, max_m, tol)
m = length(r)
resid = F - r[m].(t)
idx = locate_extrema_idx(resid, m)

figure(1)
plot(t, resid, t[idx], resid[idx], "o")
grid(true)
xlabel(L"$t$", fontsize=12)
xylims=axis()
savefig("figure13L.pdf")

supp_pt = r[m].supp_pt
max_iterations=3
tx = t[idx]
opt_r, zmin, zmax = minimax(fα, tx, supp_pt, (-1.0, 1.0), max_iterations)
resid = F - opt_r[end].(t)
figure(2)
plot(t, resid, tx, fα.(tx) - opt_r[end].(tx), "o")
grid(true)
xlabel(L"$t$", fontsize=12)
axis(xylims)
savefig("figure13R.pdf")
