using MittagLefflerFunctions
using PyPlot
using SpecialFunctions: erfcx
using Printf

α = 1/2
β = 1.0

x = range(0, 10, length=201)
y = erfcx.(x)

rvals = collect(4:2:10)
maxerr = Float64[]

figure(1)
@printf("%5s  %5s  %10s\n\n", "m", "n", "max error")
style = ("-", "--", ":", "-.")
for row = 1:length(rvals)
    r = rvals[row]
    m = r + 1
    n = r
    E = MLPade(α, β, m, n, :svd)
    err = abs.( y - E.(x) )
    semilogy(x, err, style[row], label="m, n = $m, $n")
    push!(maxerr, maximum(err))
    @printf("%5d  %5d  %10.2e\n", m, n, maxerr[row])
end
legend()
grid(true)
xlabel(L"$x$")
savefig("figure9.pdf")
