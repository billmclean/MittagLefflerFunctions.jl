using MittagLefflerFunctions
using LinearAlgebra: cond, svd, lu, norm
using PyPlot
using Printf

α = 0.2
β = 1.0
r = 8
m = r+1
n = r
C = MittagLefflerFunctions.Pade.pade_matrix(α, β, m, n)
Ctilde = [ C[:,1:r] C[:,r+2:2r+1] ]
@printf("Condition number of Ctilde = %8.2e\n", cond(Ctilde))
@printf("Condition number of C      = %8.2e\n", cond(C))
F_lu = lu(C)
@printf("Condition number of U      = %8.2e\n", cond(F_lu.U))

R_fix = MLPade(α, β, m, n, :fix)
R_svd = MLPade(α, β, m, n, :svd)
R_lu  = MLPade(α, β, m, n, :lu)

@printf("\t||p_fix-p_svd||_oo = %8.2e\n", norm(R_fix.p-R_svd.p, Inf))
@printf("\t||q_fix-q_svd||_oo = %8.2e\n", norm(R_fix.q-R_svd.q, Inf))
@printf("\t||p_lu-p_svd||_oo  = %8.2e\n", norm(R_lu.p-R_svd.p, Inf))
@printf("\t||q_lu-q_svd||_oo  = %8.2e\n", norm(R_lu.q-R_svd.q, Inf))

figure(1)
x = range(0, 5, length=201)
plot(x, R_fix.(x) - R_svd.(x), x, R_lu.(x) - R_svd.(x))
legend(("fix - svd", "lu - svd"))
grid(true)
