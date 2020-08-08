using MittagLefflerFunctions
using Polynomials
using PyPlot

α = 3/4
β = 1.0

E = MLQuad(α, β, 10, :hyperbola)

Nx = 300
Ny = 300

xmin, xmax = -6.0, 2.0
ymin, ymax = -4.0, 4.0
x = range(xmin, xmax, length=Nx)
y = range(ymin, ymax, length=Ny)

z = Complex{Float64}[ complex(x[j], y[i]) for i=1:Ny, j=1:Nx ]
modE = abs.(E.(z))

figure(1)
contourf(x, y, log10.(modE), levels=range(-2, 2, length=13))
colorbar(ticks=range(-2, 2, length=5))
xlabel(L"$\Re z$")
ylabel(L"$\Im z$")
grid(true)
savefig("figure02.pdf")
