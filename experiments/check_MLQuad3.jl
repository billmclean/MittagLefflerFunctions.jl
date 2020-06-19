using MittagLefflerFunctions
using PyPlot

α = 1.85
β = 1.0
x = range(0, 10, length=201)
E1 = MLQuad3(α, β, 6)
E2 = MLQuad3(α, β, 8)

y = mlf.(α, β, -x, 10)
y1 = E1.(x)
y2 = E2.(x)

figure(1)
#plot(x, y1-y, x, y2-y)
semilogy(x, abs.(y1-y), x, abs.(y2-y))
grid(true)
