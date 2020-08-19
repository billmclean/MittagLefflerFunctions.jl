using MittagLefflerFunctions
import MittagLefflerFunctions.Quadrature.MLQuad
using RationalApproximations
using PyPlot
using Printf

G(t) = (1-t)/(1+t)
big_α_vals = collect(1:5) ./ big(5)
α_vals = collect(1:5) ./ 5
β = big(1.0)
#big_Eα = [ MLQuad(α, β, 20, :hyperbola) for α in big_α_vals ]
Eα = [ MLQuad(α, 1.0, 14, :hyperbola) for α in α_vals ]
fα = Vector{Function}(undef, 5)
for col = 1:5
#    fα[col] = ( t -> convert(Float64, (big_Eα[col](-G(big(t))))) )
    fα[col] = ( t -> Eα[col](-G(t)) )
end
t = range(-1, 1, length=800)
max_m = 10
tol = 1e-10
max_iterations = 5

rows = 10
cols = length(α_vals)
errors = Array{Float64}(undef, rows, cols)
δ = similar(errors)
for col = 1:cols
    α = α_vals[col]
    println("α = ", α)
    F = fα[col].(t)
    r, err = aaa_v2(F, t, max_m, tol)
    M = length(r)
    for row = 1:rows
        if row ≤ M
            m = row
            resid = F - r[m].(t)
            idx = locate_extrema_idx(resid, m)
            if length(idx) ≠ 2m
                errors[row,col] = NaN
                δ[row,col] = NaN
            else
                supp_pt = r[m].supp_pt
                opt_r, zmin, zmax = minimax(fα[col], t[idx], supp_pt, 
                                    (-1.0,1.0), max_iterations) 
                errors[row,col] = zmax[end]
                δ[row,col] = zmax[end] - zmin[end]
            end
        else
            errors[row,col] = NaN
            δ[row,col] = NaN
        end
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
