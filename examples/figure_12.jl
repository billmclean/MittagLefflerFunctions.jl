using MittagLefflerFunctions
import MittagLefflerFunctions.Quadrature.MLQuad
using PyPlot

α_vals = [ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ]
β = 1.0
G(t) = (1-t)/(1+t)
t = range(-1, 1, length=400)
x = G.(t)
Eα = Vector{MLQuad}(undef, 10)

style = ("-", "--", ":", "-.", "-")

figure(1)
for j = 1:5
    α = α_vals[j]
    Eα[j] = MLQuad(α, β, 12, :hyperbola)
    plot(t, Eα[j].(-x), style[j], label="α = $α")
end
grid(true)
legend()
xlabel(L"$t$", fontsize=12)
savefig("figure12L.pdf")

figure(2)
s = range(0, 1, length=400)
t = 2 * (s.^2) .- 1
x = G.(t)
for j = 6:10
    α = α_vals[j]
    Eα[j] = MLQuad(α, β, 12, :hyperbola)
    plot(t[2:end], Eα[j].(-x[2:end]), style[j-5], label="α = $α")
end
grid(true)
legend()
xlabel(L"$t$", fontsize=12)
savefig("figure12R.pdf")
