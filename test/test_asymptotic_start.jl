import MittagLefflerFunctions.asymptotic_start
using PyPlot

β = 1.0
m = 201
α = range(0, 1, length=m)
tol = 1e-15
min_x = 35

x = Vector{Float64}(undef, m)
N = Vector{Int64}(undef, m)
for k = 1:m
    x[k], N[k] = asymptotic_start(α[k], β, tol, min_x)
end

figure(1)
subplot(2, 1, 1)
plot(α, x, "o", markersize=2)
ylabel(L"$x$")
grid(true)
subplot(2, 1, 2)
plot(α, N, "o", markersize=2)
xlabel(L"$\alpha$")
ylabel(L"$N$")
grid(true)
