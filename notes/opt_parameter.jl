using PyPlot
using Optim

A(ϕ) = acosh( 2π*ϕ / ( (4ϕ-π) * π * sin(ϕ) ) )
B(ϕ) = ( π^2 - 2π*ϕ ) / A(ϕ)

res = optimize(ϕ -> -B(ϕ), π/4, π/2)
ϕ_opt = res.minimizer

ϕ = range(π/4, π/2, length=201)
figure(1)
plot(ϕ, B.(ϕ), [ϕ_opt], [B(ϕ_opt)], "o")
grid(true)
xlabel("ϕ")
savefig("Bplot.pdf")

println("Optimal values:")
println("\tϕ = ", ϕ_opt)
println("\th = $(A(ϕ_opt)) / N")
println("\tμ = $((4π*ϕ_opt-π^2)/A(ϕ_opt)) N / t")
println("\tB(ϕ) = ", B(ϕ_opt))
println("\tDecay factor = $(exp(B(ϕ_opt)))^(-N)")

