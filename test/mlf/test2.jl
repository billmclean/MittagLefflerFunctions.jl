import MittagLefflerFunctions.MLF.mlf2
import MittagLefflerFunctions.MLF.radius

α = 0.5
β = 1.0
z = 8.7 * exp(Complex(0.0, 3π/4))
p = 10
@assert (   (  abs(angle(z)) > π*α )
         && ( abs(abs(angle(z)) - (π*α)) > 10.0^(-p) ) )

r0 = radius(α, β, z, p)
val = mlf2(α, β, z, p, r0)
ref = erfcx(-z)

@test abs(val-ref) ≤ abs(ref)
