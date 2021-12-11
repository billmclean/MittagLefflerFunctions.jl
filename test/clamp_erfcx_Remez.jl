import PolynomialApproximations: Remez
using SpecialFunctions: erfcx
using PyPlot

T = Float64

function g(y::T) where T <: AbstractFloat
    if -one(T) ≤ y ≤ -one(T) + eps(T)
        return one(T)
    elseif one(T) - eps(T) ≤ y ≤ one(T)
        return zero(T)
    else
        x = (1+y) / (1-y)
        return erfcx(x)
    end
end

function erfcx_approx(x)
    if x < one(T)
        y = (x-1) / (x+1)
    else
        rx = 1 / x
        y = (1-rx) / (1+rx)
    end
    return g_approx(y)
end

n = 15
iterations = 6
clamp = :both
ddp, pt, zmax, zmin = Remez.minimax(T, g, n, iterations, clamp)

g_approx(y) = Remez.Newton_poly(ddp, pt, y)

figure(1)
y = range(-one(T), one(T), length=201)

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
x = Float64[ 10 * (n/400)^2 for n = 0:400 ]
plot(x, erfcx_approx.(x) - erfcx.(x))
grid(true)
xlabel(L"$x$")
ylabel("error in erfcx")
