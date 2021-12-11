import MittagLefflerFunctions.H
import MittagLefflerFunctions.Fplus
using PyPlot

α = 0.7
β = 0.95
tol = 1e-10

x = range(0.5, 1.5, length=201)
z = complex.(x)
errH(z) = (z^α-1) - (z-1)*(α+(z-1)*H(α, z, tol))

y = real(errH.(z))
figure(1)
plot(x, y)
grid(true)

G(z) = z^(α-β) * (z-1) / (z^α-1)
errG(z) = G(z) - z^(α-β) / (α+(z-1)*H(α,z,tol))
y = real(errG.(z))
figure(1)
figure(2)
plot(x, y)
grid(true)

y = real(Fplus.(α, β, z, tol))
figure(3)
lo = 1 - 0.1
hi = 1 + 0.1
plot(x, y)
v = axis()
plot([lo, lo, NaN, hi, hi], [v[3], v[4], NaN, v[3], v[4]], ":")
grid(true)
