import MittagLefflerFunctions.MLF.mlf1

α = 1/2
β = 1.0
z = Complex(-0.5, 0.4)
@assert ( α < 1 ) && ( abs(z) ≤ 1 )

val = mlf1(α, β, z, 10)
ref = erfcx(-z)

@test abs(val-ref) ≤ 1e-10

α = 1.5
β = 1.0
z = Complex(5.5, 0.0)
@assert ( 1 ≤ α < 2 ) && ( abs(z) ≤ floor(20/(2.1-α)^(5.5-2α)) )

val = mlf1(α, β, z, 10)
Eq = MLQuad2(α, β, 10)
ref = Eq(real(z))

@test abs(val-ref) ≤ 1e-10 * abs(ref)

α = 4.0
β = 1.0
z = Complex(15.0, 17.3)
@assert ( α ≥ 2 ) && ( abs(z) ≤ 50 ) 

val = mlf1(α, β, z, 10)
ref = ( cos(z^(1/4)) + cosh(z^(1/4)) ) / 2

@test abs(val-ref) ≤ 1e-10 * abs(ref)

