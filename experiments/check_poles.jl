import MittagLefflerFunctions.Quadrature.H2
import MittagLefflerFunctions.Quadrature.H3
using PyPlot

α = 1.7
β = 1.2
sep = 0.02
x = 2.0
M = 200
t = range(-0.2, 0.2, length=M)
θ = π/4

w = x^(1/α) .+ t * exp(Complex(0.0, θ))
y = H2.(α, β, w, x, sep) 

figure(1)
H0 = ( 1 + α - 2β ) / ( 2α * x^(β/α) )
plot(t, real.(y), t, imag.(y), [0], [H0], "o")
grid(true)

γ₊ = x^(1/α) * exp(Complex(0.0, π/α))
γ₋ = x^(1/α) * exp(Complex(0.0, -π/α))
w = γ₊ .+ t * exp(Complex(0.0, θ))

y1 = zeros(Complex{Float64}, M)
y2 = zeros(Complex{Float64}, M)
y3 = zeros(Complex{Float64}, M)
for m = 1:M
    H₊, H₋ = H3(α, β, w[m], x, sep)
    y1[m] = H₊
    y2[m] = H₋
    y3[m] = ( ( w[m] - γ₊ ) * ( w[m] - γ₋ ) 
             / ( ( w[m]^α + x ) * ( 2w[m] - γ₊ - γ₋ ) ) )
end

figure(2)
H0₊ = ( ( ( 1 + α - 2β ) * ( γ₊ - γ₋ ) - 2γ₊ )
      / ( 2α * γ₊^β * ( γ₊ - γ₋ ) ) )
plot(t, real.(y1), t, imag.(y1),
     [0], [real(H0₊)], "o", [0], [imag(H0₊)], "o")
grid(true)

figure(3)
H0₋ = ( sinpi((1-β)/α) / ( α * sinpi(1/α) ) ) / x^(β/α)
plot(t, real.(y2), t, imag.(y2), [0], [H0₋], "o")
grid(true)

figure(4)
ϕ0 = γ₊^(1-α) / α
plot(t, real.(y3), t, imag.(y3), 
     [0], [ real(ϕ0) ], "o", [0], [ imag(ϕ0) ], "o")
grid(true)
