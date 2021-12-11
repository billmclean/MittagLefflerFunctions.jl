using SpecialFunctions: loggamma
using Printf: @printf

T = BigFloat
case = 1
if case == 1
    α_vals = T[ n/10 for n = 1:2:9 ]
    x_vals = T[ n*5 for n = 1:8 ]
    β = 13 / parse(T, "17")
elseif case == 2
    α_vals = T[ (994+n)/1000 for n = 1:5 ]
    x_vals = T[ n*5 for n = 5:10 ]
    β = one(T)
end
@printf("\nβ = %10.6f\n", β)

tol = parse(T, "10")^(-16)
include("assess.jl")

@printf("\nMinimum N such that Nth term is less than %0.2e\n\n", tol)
@printf("  x    ")
for α in α_vals
    @printf("   α = %5.3f  ", α)
end
@printf("\n\n")
for x in x_vals
    @printf("%5.2f | ", x)
    for α in α_vals
        Nmin, term = assess(α, β, x, tol)
        if term < tol
            @printf("%3d %8.2e |", Nmin, term)
        else
            @printf("%3s %8s |", "**", "**")
        end
    end
    @printf("\n")
end
@printf("\n\nOptimal truncation at N = x^(1/α):\n")
@printf("\n  x    ")
for α in α_vals
    @printf("   α = %5.3f  ", α)
end
@printf("\n\n")
for x in x_vals
    @printf("%5.2f | ", x)
    for α in α_vals
        Nopt = ceil(Integer, x^(1/α))
        if Nopt > 10_000_000_000
            @printf(" %11s |", "**")
        else
            @printf(" %11d |", Nopt)
        end
    end
    @printf("\n")
end
