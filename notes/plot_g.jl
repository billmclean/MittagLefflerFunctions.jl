using PyPlot
using SpecialFunctions
using MittagLefflerFunctions
using OffsetArrays

function g(y)
    x = (1+y) / (1-y)
    return ( erfcx(x) - one(y) ) / x
end

function Eneg(a, x)
    y = (x-1) / (x+1)
    return 1 + x * chebyshev_sum(a, y)
end

nmax = 20
a = OffsetArray{Float64}(undef, 0:nmax)
M = 2nmax
chebyshev_coefs!(a, g, M)

figure(1)
y = range(-1, 1, length=201)
plot(y, g.(y))
grid(true)

figure(2)
approx_g = Float64[ chebyshev_sum(a, yval) for yval in y ]
plot(y, approx_g - g.(y))
grid(true)

figure(3)
x = range(0, 5, length=201)
E = Float64[ Eneg(a, xval) for xval in x ]
plot(x, E - erfcx.(x))
grid(true)
 
