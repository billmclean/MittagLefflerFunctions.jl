import Remez
using SpecialFunctions: erfcx
using PyPlot

T = Float64

function g(y::T) where T <: AbstractFloat
    two = parse(T, "2")
    rootpi = sqrt(convert(T, π))
    if -one(T) ≤ y ≤ -one(T) + eps(T)
        return ( 1 / two - 1 / rootpi )
    elseif one(T) - eps(T) ≤ y ≤ one(T)
        return ( -1 + 1/rootpi ) / 2
    else
        x = (1+y) / (1-y)
        return ( 2*erfcx(x) / (1-y) - 1 ) / (1+y)
    end
end

function erfcx_approx(x)
    y = (x-1) / (x+1)
    return (1-y)*(1+(1+y)*g_approx(y)) / 2
end

n = 15
iterations = 6
pt, ddp, E = Remez.minimax(T, g, n, iterations)

g_approx(y) = Remez.Newton_poly(ddp[0:n], pt[0:n], y)

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
x = range(0, 10, length=201)
plot(x, erfcx_approx.(x) - erfcx.(x))
grid(true)
xlabel(L"$x$")
ylabel("error in erfcx")
