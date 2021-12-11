using MittagLefflerFunctions: ML
using PyPlot
using SpecialFunctions: erfcx

α = 0.5
β = 1.0
tol = 1e-10

x = range(-3, 3, length=501)
y1 = ML.(α, β, x, tol)
y2 = erfcx.(-x)

figure(1)
semilogy(x, y1, x, y2, "--")
grid(true)

figure(2)
plot(x, abs.(y1-y2)./y2)
grid(true)
