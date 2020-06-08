module Quadrature

using OffsetArrays

export MLQuad1, MLQuad2, MLQuad3

const ϕ_opt = "1.17210"

coef_h(ϕ) = acosh( 2ϕ / ((4ϕ-π)*sin(ϕ)) )
coef_μ(ϕ) =  π*(4ϕ-π) / coef_h(ϕ)

"""
`QSum` holds the coefficients for a quadrature sum approximating the 
contour integral of `e^w F(w) / 2πi`:

    A * ( C[0] Re( F(w[0]) ) + ⋯ + C[N] * Re( F(w[N]) ) )
"""
struct QSum{T<:AbstractFloat}
    w::OffsetArray{Complex{T}}
    C::OffsetArray{Complex{T}}
    A::T
end

"""
`MLQuad1` holds data for computing the Mittag--Leffler function `Eαβ(-x)`
when `x≥0` and `0<α<1`.
"""
struct MLQuad1{T<:AbstractFloat}
    α::T
    β::T
    qs::QSum{T}
end

"""
`MLQuad2` holds data for computing the Mittag--Leffler function `Eαβ(x)`
when `x≥0` and `0<α<2`.
"""
struct MLQuad2{T<:AbstractFloat}
    α::T
    β::T
    qs::QSum{T}
end

"""
`MLQuad3` holds data for computing the Mittag--Leffler function `Eαβ(-x)`
when `x≥0` and `1<α<2`.
"""
struct MLQuad2{T<:AbstractFloat}
    α::T
    β::T
    qs::QSum{T}
end

function QSum(::Type{T}, N) where T <: AbstractFloat
    ϕ = parse(T, ϕ_opt)
    a, b = coef_h(ϕ), coef_μ(ϕ)
    h = a / N
    μ = b * N
    w = OffsetArray{Complex{T}}(undef, 0:N)
    C = OffsetArray{Complex{T}}(undef, 0:N)
    w0 = μ * ( 1 - sin(ϕ) )
    w[0] = Complex(w0, 0)
    C[0] = Complex(exp(w0) * cos(ϕ) / 2, 0)
    for n = 1:N
        un = n * h
        iunmϕ = Complex(-ϕ, un)
        w[n] = μ * ( 1 + sin(iunmϕ) )
        C[n] = exp(w[n]) * cos(iunmϕ)
    end
    A = a * b / π
    return QSum(w, C, A)
end

function MLQuad1(α::T, β::T, N::Integer) where T <: AbstractFloat
    if !(0≤α≤1)
        throw(DomainError(α, "α must lie between 0 and 1"))
    end
    qs = QSum(T, N)
    return MLQuad1(α, β, qs)
end

function (E::MLQuad1{T})(x::T) where T <: AbstractFloat
    α, β, qs = E.α, E.β, E.qs
    w, C, A = qs.w, qs.C, qs.A
    if x < 0
        throw(DomainError(x, "argument must be greater than or equal to zero"))
    else
        s = zero(T)
        for n in axes(w, 1)
            s += real( C[n] * w[n]^(α-β) / ( x + w[n]^α ) )
        end
        return A * s
    end
end

function MLQuad2(α::T, β::T, N::Integer) where T <: AbstractFloat
    if !(0≤α≤2)
        throw(DomainError(α, "α must lie between 0 and 2"))
    end
    qs = QSum(T, N)
    return MLQuad2(α, β, qs)
end

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

function (E::MLQuad2{T})(x::T) where T <: AbstractFloat
    α, β, qs = E.α, E.β, E.qs
    w, C, A = qs.w, qs.C, qs.A
    sep = 1 / parse(T, "4")
    if x < 0
        throw(DomainError(x, "argument must be greater than or equal to zero"))
    else
        s = zero(T)
        for n in axes(w, 1)
            s += real( C[n] * H(α, β, w[n], x, sep) )
        end
        return ( x^((1-β)/α) * exp(x^(1/α)) / α ) + A * s
    end
end

function MLQuad3(α::T, β::T, N::Integer) where T <: AbstractFloat
    if !(1≤α≤2)
        throw(DomainError(α, "α must lie between 1 and 2"))
    end
    qs = QSum(T, N)
    return MLQuad3(α, β, qs)
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

function (E::MLQuad3{T})(x::T) where T <: AbstractFloat
    α, β, qs = E.α, E.β, E.qs
    w, C, A = qs.w, qs.C, qs.A
    sep = 1 / parse(T, "4")
    if x < 0
        throw(DomainError(x, "argument must be greater than or equal to zero"))
    else
        s = zero(T)
        for n in axes(w, 1)
            s += real( C[n] * ( Hplus(α, β, w[n], x, sep) 
                             + Hminus(α, β, w[n], x, sep) ) ) / 2
        end
        sum_residues = ( ( x^((1-β)/α) / α ) * exp(x^(1/α)*cos(π/α))
                 * cos( π*(1-β)/α + x^(1/α)*sin(π/α) ) )
        return sum_residues + A * s
    end
end

end # module
