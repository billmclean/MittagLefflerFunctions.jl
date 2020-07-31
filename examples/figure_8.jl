using MittagLefflerFunctions
using SpecialFunctions: erfcx
using PyPlot

α = 1/2
β = 1.0
r = 5

x = range(0, 10, length=201)
y = erfcx.(x)

figure(1)
m_lo = 4
style = ("-", "--", ":", "-.", "-")
for m = m_lo:2r+2-m_lo
    n = 2r + 1 - m
    E = MLPade(α, β, m, n, :svd)
    err = E.(x) - y
    plot(x, err, style[m-m_lo+1], label="m = $m, n = $n")
end
xlabel(L"$x$")
legend()
grid(true)

savefig("figure8.pdf")
