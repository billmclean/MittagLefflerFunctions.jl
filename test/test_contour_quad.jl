import MittagLefflerFunctions.MLneg_integral
using Printf
using PyPlot

T = BigFloat
lenA = 20
A = Vector{T}(undef, lenA)
α = parse(T, "0.6")
β = one(T)
t = parse(T, "0.1")
λ = one(T)

for N = 1:lenA
    A[N] = MLneg_integral(α, β, λ, t, N)
    if N ≥ 2
        @printf("%5d  %14.4e\n", N, A[N]-A[N-1])
    end
end

figure(1)
semilogy(collect(1:lenA-5), abs.(A[1:end-5].-A[end]), "o")
grid(true)
