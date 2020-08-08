using MittagLefflerFunctions
using Polynomials
using PyPlot
using SpecialFunctions: erfcx

α = 1/2
β = 1.0
N = 5
r = 2N + 1
m = r + 1
n = r

#Eref = MLQuad(α, β, 2N, :hyperbola)
Eref(z) = erfcx(-z)

EP = MLPade(α, β, m, n, :svd)
p, q = Polynomial(EP.p), Polynomial(EP.q)
χ = roots(q)
EPade(z) = p(-z) / q(-z)

EQ = MLQuad(α, β, N, :hyperbola)
w, C, A = EQ.qs.w, EQ.qs.C, EQ.qs.A
wα = w.^α
P = [ conj.(wα[N:-1:1]); wα[0:N] ]
R1 = -A * C .* w.^(α-β)
R = [ conj.(R1[N:-1:1]); R1[0:N] ]

function mlf_quad1(z)
    s = 0.0
    for j = 1:length(P)
        s += R[j] / ( z - P[j] )
    end
    return s
end

figure(1)
plot(real.(-χ), imag.(-χ), "o",
     real.(P), imag.(P), "o")
legend((L"$-\chi_j$", L"$P_j$"))
grid(true)

Nx = 301
Ny = 301
xmin, xmax = -5.0, 3.0
ymin, ymax = -4.0, 4.0
x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)

z = Complex{Float64}[ complex(x[j], y[i]) for i=1:Ny, j=1:Nx ]

err_quad = Float64[ abs(Eref(z[i,j])-mlf_quad1(z[i,j])) for i = 1:Ny, j=1:Nx ]
err_pade = Float64[ abs(Eref(z[i,j])-EPade(z[i,j])) for i = 1:Ny, j=1:Nx ]

figure(2)
contour(x, y, log10.(err_quad), levels=collect(-9:4))
plot(real.(P), imag.(P), "ro")
colorbar()
xlabel(L"$\Re z$")
ylabel(L"$\Im z$")
grid(true)
savefig("figure11L.pdf")

figure(3)
contour(x, y, log10.(err_pade), levels=collect(-14:4))
plot(real.(-χ), imag.(-χ), "ro")
colorbar()
xlabel(L"$\Re z$")
ylabel(L"$\Im z$")
grid(true)
savefig("figure11R.pdf")
                        
