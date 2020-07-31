import MittagLefflerFunctions.Quadrature.MLQuad
using SpecialFunctions: erfcx
using PyPlot

xmin, xmax = -0.4, 0.4
ymin, ymax = -2.5, 2.5
Nx = 201
Ny = 201
x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)

z = Complex{Float64}[ complex(x[j], y[i]) for i=1:Ny, j=1:Nx ]

E_par = MLQuad(0.5, 1.0, 10, :parabola, 0.2)
E_hyp = MLQuad(0.5, 1.0, 10, :hyperbola, 0.2)

err_par = Float64[ abs(E_par(z[i,j])-erfcx(-z[i,j])) for i=1:Ny, j=1:Nx ]
err_hyp = Float64[ abs(E_hyp(z[i,j])-erfcx(-z[i,j])) for i=1:Ny, j=1:Nx ]

figure(1)
#mesh(x, y, err_par)
contour(x, y, err_par, 12)
xlabel(L"$\Re z$")
ylabel(L"$\Im z$")
grid(true)
colorbar()
savefig("figure6L.pdf")

figure(2)
contour(x, y, err_hyp, 12)
xlabel(L"$\Re z$")
ylabel(L"$\Im z$")
grid(true)
colorbar()
savefig("figure6R.pdf")
