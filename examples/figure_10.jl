using MittagLefflerFunctions
using PyPlot
using Polynomials

x = range(0, 10, length=201)
α_vals = [0.2, 0.4, 0.6, 0.8, 1.0]

figure(1)
r = 5
m = r+1
n = r
style = ("-", "--", ":", "-.", "-")
for row = 1:length(α_vals)
    α = α_vals[row]
    Eα_approx = MLPade(α, 1.0, m, n, :svd)
    Eα_ref = MLQuad(α, 1.0, 14, :hyperbola)
    labelstring = latexstring("\$\\alpha=$α\$")
    semilogy(x, abs.(Eα_ref.(-x) - Eα_approx.(x)), 
             style[row], label=labelstring)
end
xlabel(L"$x$")
legend()
grid(true)

savefig("figure10.pdf")

α = 1.85
β = 1.0
println("Setting α = $α, β = $β.")
for (m, n) in [(8,7), (10,9)]
    println("\nRoots of q when m = $m, n = $n:")
    E = MLPade(α, β, m, n, :svd)
    q = Polynomial(E.q)
    display(roots(q))
end
