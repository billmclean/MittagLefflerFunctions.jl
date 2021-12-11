using SpecialFunctions: loggamma
using PyPlot

include("assess.jl")

α_vals = Float64[ n/50 for n = 10:2:50 ]
β = 1.0
tol = 1e-16
x = range(1, 50, length=50)
Nmin = Array{Float64}(undef, 50)

figure(1)

for α in α_vals
    for k = 1:50
        Nmin[k], term = assess(α, β, x[k], tol)
        if Nmin[k] < 0
            Nmin[k] = NaN
        end
    end
    plot(x, Nmin, label="α = $α")
end
title("tol = $tol")
#legend()
xlabel("x")
ylabel("N")
grid(true)
   
