using MittagLefflerFunctions
import RationalApproximations: locate_extrema_idx, aaa_v2, minimax
using PyPlot

α = 0.75
β = 1.0
N = 4
Equad = MLQuad(α, β, N, :hyperbola)
E_ref = MLQuad(α, β, 14, :hyperbola)

function f(t)
    x = (1-t) / (1+t)
    return E_ref(-x)
end

t = collect(range(-1, 1, length=1000))
max_m = 2N+1
tol = 1e-14
F = f.(t)
r, err = aaa_v2(F, t, max_m, tol)
#m = length(r)
m = 2N
resid = F - r[m].(t)
idx = locate_extrema_idx(resid, m)

figure(1)
plot(t, resid, t[idx], resid[idx], "o")
xlabel(L"t")
grid(true)
title("Error in AAA approximation with m = $m")

max_iterations = 5
clamp = :no
opt_r, zmin, zmax = minimax(f, t[idx], r[m].supp_pt, (-1.0,1.0), 
                             max_iterations, clamp)
opt_resid = F - opt_r[end].(t)

figure(2)
x = (1 .- t)./(1 .+ t)
Equad_resid = Equad.(-x) - F
plot(t[1:end-1], Equad_resid[1:end-1])
grid(true)
xlabel(L"$t$")
title("Error in quadrature approximation with N = $N")

figure(3)
plot(t, opt_resid, "-")
legend()
xlabel(L"t")
title("Error in best rational approximation with m = $m")
grid(true)

figure(4)
loglog(x[1:end-1], abs.(Equad_resid[1:end-1]), "-", label="N = $N")
loglog(x[1:end-1], abs.(opt_resid[1:end-1]), "--", label="m = $m")
#axis([-1.1, 1.1, 5e-11, 5e-5])
legend()
xlabel(L"t")
grid(true)
savefig("figure14.pdf")

