using PyPlot
using Optim

A(ϕ) = acosh( 2ϕ / ( (4ϕ-π) * sin(ϕ) ) )
B(ϕ) = ( π^2 - 2π*ϕ ) / A(ϕ)

res = optimize(ϕ -> -B(ϕ), π/4, π/2)
ϕ_opt = res.minimizer

δ = 1e-3
ϕ = range(π/4+δ, π/2-δ, length=201)
figure(1)
plot(ϕ, B.(ϕ), [ϕ_opt], [B(ϕ_opt)], "o")
xticks((π/4, 3π/8, π/2),(L"$\pi/4$", L"$3π/8$", L"$\pi/2$"),fontsize=12)
grid(true)
savefig("Bplot.pdf")

println("Optimal values:")
println("\tϕ = ", ϕ_opt)
println("\th = $(A(ϕ_opt)) / N")
println("\tμ = $((4π*ϕ_opt-π^2)/A(ϕ_opt)) N")
println("\tB(ϕ) = ", B(ϕ_opt))
println("\tDecay factor = $(exp(B(ϕ_opt)))^(-N)")

