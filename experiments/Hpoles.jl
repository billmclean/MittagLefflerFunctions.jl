using PyPlot
using SpecialFunctions: erfcx

"""
    ψ1(α, ϵ, N=100)

Approximates `(1+ϵ)^α - 1` by its Taylor expansion about `ϵ=0` using up to 
`N` terms.
"""
function ψ1(α::T, ϵ::Complex{T}, N=100) where T <: AbstractFloat
    s = Complex(α, zero(T))
    term = s
    for n = 1:N
        term *= -( (n-α)/(n+1) ) * ϵ
        s += term
        if abs(term) < eps(T)
            break
        end
    end
    if abs(term) > 10*eps(T)
        error("ψ1 expansion did not converge after $N terms")
    end
    return s
end

"""
    ψ2(α, ϵ, N=100)

Approximates `(1+ϵ)^α - (1+αϵ)` by its Taylor expansion about `ϵ=0` using up 
to `N` terms.
"""
function ψ2(α::T, ϵ::Complex{T}, N=100) where T <: AbstractFloat
    s = Complex(-α * (1-α) / 2, zero(T))
    term = s
    for n = 1:N
        term *= -( (n+1-α)/(n+2) ) * ϵ
        s += term
        if abs(term) < eps(T)
            break
        end
    end
    if abs(term) > 10*eps(T)
        error("ψ1 expansion did not converge after $N terms")
    end
    return s
end

function H(α::T, β::T, w::Complex{T}, x::T, sep::T) where T <: AbstractFloat
    xa = x^(1/α)
    ϵ = ( w - xa ) / xa
    if abs(ϵ) > sep
        return w^(α-β) / (w^α-x) - 1 / (α*ϵ*x^(β/α))
    else
        return ( ψ1(α-β,ϵ) - (1+ϵ)^(α-β) * ψ2(α,ϵ)/ψ1(α,ϵ) ) / (α*x^(β/α))
    end
end

    
a(ϕ) = acosh(2ϕ/((4ϕ-π)*sin(ϕ)))
b(ϕ) = ( (4π*ϕ-π^2) / a(ϕ) ) 

function Qpos(α::T, β::T, x::T, N::Integer, sep::T) where T <: AbstractFloat
    ϕ = parse(T, "1.17210")
    h = a(ϕ) / N
    μ = b(ϕ) * N
    w0r = μ * ( 1 - sin(ϕ) )
    w0 = Complex(w0r, zero(T))
    s = ( exp(w0r) * cos(ϕ) / 2 ) * real(H(α, β, w0, x, sep)) 
    for n = 1:N
        un = n * h
        iunmϕ = Complex(-ϕ, un)
        wn = μ * ( 1 + sin(iunmϕ) )
        s += real( exp(wn) * H(α, β, wn, x, sep) * cos(iunmϕ) )
    end
    return ( x^((1-β)/α) * exp(x^(1/α)) / α ) + ( a(ϕ) * b(ϕ) / π ) * s
end

function Hplus(α::T, β::T, w::Complex{T}, x::T, sep::T) where T <: AbstractFloat
    γp = x^(1/α) * exp(complex(zero(T), π/α))
    ϵp = ( w - γp ) / γp
    if abs(ϵp) > sep
        return w^(α-β) / ( w^α + x ) - γp^(1-β) / ( α * ϵp * γp^β )
    else
        return ( ψ1(α-β,ϵp) - (1+ϵp)^(α-β)*ψ2(α,ϵp) / ψ1(α,ϵp) ) / (α*γp^β)
    end
end

function Hminus(α::T, β::T, w::Complex{T}, x::T, sep::T
               ) where T <: AbstractFloat
    γm= x^(1/α) * exp(complex(zero(T), -π/α))
    ϵm = ( w - γm) / γm
    if abs(ϵm) > sep
        return w^(α-β) / ( w^α + x ) - γm^(1-β) / ( α * ϵm * γm^β )
    else
        return ( ψ1(α-β,ϵm) - (1+ϵm)^(α-β)*ψ2(α,ϵm) / ψ1(α,ϵm) ) / (α*γm^β)
    end
end

function Qneg(α::T, β::T, x::T, N::Integer, sep::T) where T <: AbstractFloat
    @assert 1 < α < 2
    ϕ = parse(T, "1.17210")
    h = a(ϕ) / N
    μ = b(ϕ) * N
    w0r = μ * ( 1 - sin(ϕ) )
    w0 = Complex(w0r, zero(T))
    s = ( exp(w0r) * cos(ϕ) / 2 ) * real(Hplus(α, β, w0, x, sep)) 
    for n = 1:N
        un = n * h
        iunmϕ = Complex(-ϕ, un)
        wn = μ * ( 1 + sin(iunmϕ) )
        s += real( exp(wn) * ( Hplus(α, β, wn, x, sep) 
                             + Hminus(α, β, wn, x, sep) ) * cos(iunmϕ) ) / 2
    end
    sum_residues = ( ( x^((1-β)/α) / α ) * exp(x^(1/α)*cos(π/α)) 
                     *cos( π*(1-β)/α + x^(1/α)*sin(π/α) ) ) 
    return sum_residues + ( a(ϕ)*b(ϕ)/π ) * s
end

α = 0.5
β = 1.0
sep = 0.2
N = 10

figure(1)
x = 2.0
t = range(-2.0, 2.0, length=201)
θ = π/4
w = x^(1/α) .+ t * exp(Complex(0.0, θ))
Hvals = H.(α, β, w, x, sep) 
plot(t, real.(Hvals), t, imag.(Hvals))
grid(true)

figure(2)
x = range(0.0, 1.0, length=201)
plot(x, Qpos.(α, β, x, N, sep) - erfcx.(-x))
grid(true)

figure(3)
α = 3/2
x = 2.0
γplus = x^(1/α) * exp(Complex(0.0,π/α))
t = range(-2.0, 2.0, length=201)
θ = π/4
w = x^(1/α) .+ t * exp(Complex(0.0, θ))
Hvals = Hplus.(α, β, w, x, sep) 
plot(t, real.(Hvals), t, imag.(Hvals))
grid(true)
     
figure(4)
x = range(0.0, 5.0, length=201)
plot(x, Qneg.(α, β, x, N, sep))
grid(true)
