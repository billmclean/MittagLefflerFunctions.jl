using MittagLefflerFunctions
using PyPlot

f(x) = exp(x)

nmax = 10
M = nmax + 2
a = chebyshev_coefs(Float64, f, nmax, M)

x = range(-1, 1, length=201)
y = Float64[ chebyshev_sum(a, xval) for xval in x ]

figure(1)
plot(x, y-f.(x))
grid(true)
