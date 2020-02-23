using MittagLefflerFunctions

f(x) = exp(x)

nmax = 10
a = OffsetArray{Float64}(undef, 0:nmax)
M = nmax + 2
chebyshev_coefs!(a, f, M)

x = range(-1, 1, length=201)
y = Float64[ chebyshev_sum(a, xval) for xval in x ]

figure(1)
plot(x, y-f.(x))
grid(true)
