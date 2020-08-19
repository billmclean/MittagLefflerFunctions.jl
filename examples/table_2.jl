using MittagLefflerFunctions
import MittagLefflerFunctions.Quadrature.MLQuad
using RationalApproximations
using PyPlot
using Printf

G(t) = (1-t)/(1+t)
α_vals = collect(1:5) ./ 5
Eα = [ MLQuad(α, 1.0, 14, :hyperbola) for α in α_vals ]
fα = Vector{Function}(undef, 5)
for col = 1:5
    fα[col] = ( t -> Eα[col](-G(t)) )
end
t = range(-1, 1, length=800)
max_m = 10
tol = 1e-15

rows = 10
cols = length(α_vals)
errors = fill(NaN, rows, cols)
for col = 1:cols
    α = α_vals[col]
    println("α = ", α)
    F = fα[col].(t)
    r, err = aaa_v2(F, t, max_m, tol)
    display(err)
    M = length(r)
    for row = 1:M
        errors[row,col] = err[row]
    end
end

@printf("%5s ", "m")
for col = 1:cols
    @printf("&    α = %3.1f", α_vals[col])
end
@printf("\n\n")

for row = 1:rows
    @printf("%5d ", row)
    for col = 1:cols
        @printf("& %10.2e", errors[row,col])
    end
    @printf("\\\\\n")
end

