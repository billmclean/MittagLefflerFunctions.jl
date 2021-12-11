using MittagLefflerFunctions
using SpecialFunctions: erfcx
using PyPlot

T = BigFloat

function g(y)
    x = (1+y) / (1-y)
    return ( 2*erfcx(x) / (1-y) - 1 ) / (1+y)
end

function erfcx_approx(x)
    y = (x-1) / (x+1)
    return (1-y)*(1+(1+y)*g_approx(y)) / 2
end

nmax = 30
M = nmax + 6
a = chebyshev_coefs(T, g, nmax, M)

y = range(-one(T), one(T), length=201)
g_approx(y) = chebyshev_sum(a, y)

figure(1)
plot(y, g.(y))
grid(true)
xlabel(L"$y$")
ylabel(L"$g(y)$")

figure(2)
plot(y, g.(y)-g_approx.(y))
grid(true)
xlabel(L"$y$")
ylabel(L"error in $g$")

figure(3)
semilogy(abs.(a))
grid(true)
xlabel(L"$n$")
ylabel(L"$|a_n|$")

figure(4)
x = range(zero(T), big(10.0), length=201)
plot(x, erfcx_approx.(x) - erfcx.(x))
grid(true)
xlabel(L"$x$")
ylabel("error in erfcx")
