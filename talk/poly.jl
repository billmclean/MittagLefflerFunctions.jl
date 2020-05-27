import PolynomialApproximations.Remez.equi_oscillation!
import PolynomialApproximations.NewtonPoly
using PolynomialApproximations: minimax
using SpecialFunctions: erfcx
using OffsetArrays
using Optim: optimize
using PyPlot

function g(y)
    x = ( 1 + y ) / ( 1 - y )
    return erfcx(x)
end

r = 4
pt_ = [ cospi((r+1-j)/(r+1)) for j = 0:r+1 ]
pt = OffsetArray(pt_, 0:r+1)
ddg = OffsetArray{Float64}(undef, 0:r)
ddp = OffsetArray{Float64}(undef, 0:r)

E = equi_oscillation!(ddp, pt, ddg, g, :both)
p = NewtonPoly(ddp[0:r], pt[0:r])

y = range(-1, 1, length=201)

figure(1)
err = g.(y) - p.(y)
plot(y, g.(y) - p.(y), pt, g.(pt) - p.(pt), "o")
title(latexstring("\$g_{1/2}(y)-p(y)\$ when \$r=$r\$"))
xlabel(L"$y$")
grid(true)
savefig("remez1.pdf")

new_pt = OffsetArray{Float64}(undef, 0:r+1)
z = Vector{Float64}(undef, r)
pow = ( E > 0 ) ? 1.0 : -1.0
new_pt[0] = pt[0]
for j = 1:r
    global pow
    res = optimize(pt[j-1], pt[j+1]) do y
        z = g(y) - p(y)
        return pow * z
    end
    new_pt[j] = res.minimizer
    z[j] = abs( g(new_pt[j]) - p(new_pt[j] ) )
    pow = -pow
end
new_pt[r+1] = pt[r+1]
zmax = maximum(z)
zmin = minimum(z)


figure(2)
err = g.(y) - p.(y)
plot(y, g.(y) - p.(y), new_pt, g.(new_pt) - p.(new_pt), "o")
plot([-1, 1], 
     [ zmin  zmax  -zmin  -zmax
       zmin  zmax  -zmin  -zmax ], "r", linewidth=1/2)
title(latexstring("\$g_{1/2}(y)-p(y)\$ when \$r=$r\$"))
xlabel(L"$y$")
grid(true)
savefig("remez2.pdf")

max_iterations = 10

rmin = 3
rmax = 30
rvals = OffsetArray(collect(rmin:rmax), rmin:rmax)
Evals = OffsetArray{Float64}(undef, rmin:rmax)
for r in rvals
    ddp, pt, zmax, zmin = minimax(Float64, g, r, max_iterations)
    Evals[r] = zmax[end]
end

rlo = 10
rhi = 25
A = [ ones(rhi-rlo+1) rvals[rlo:rhi] ]
b = log.(Evals[rlo:rhi])
v = A \ b
C = exp(v[1])
ρ = exp(v[2])
figure(3)
semilogy(rvals, Evals, "o", [rlo, rhi], C * [ ρ^rlo, ρ^rhi ] )
grid(true)
xlabel(L"$r$")
savefig("remez3.pdf")

println("Minimax error ≈ C ρ^r where")
println("\tC = ", C)
println("\tρ = ", ρ)

