module Quadrature

using OffsetArrays

using ..MittagLefflerFunctions: rΓ

export MLQuad, MLQuad1, MLQuad2, MLQuad3

const ϕ_opt = "1.17210"

coef_h(ϕ) = acosh( 2ϕ / ((4ϕ-π)*sin(ϕ)) )
coef_μ(ϕ) =  π*(4ϕ-π) / coef_h(ϕ)

"""
`QSum` holds the coefficients for a quadrature sum approximating the 
contour integral of `e^w F(w) / 2πi`:

    A * ( C[0] Re( F(w[0]) ) + ⋯  + C[N] * Re( F(w[N]) ) )
"""
struct QSum{T<:AbstractFloat}
    w::OffsetArray{Complex{T}}
    C::OffsetArray{Complex{T}}
    A::T
end

"""
`MLQuad` holds data for computing the Mittag--Leffler function `Eαβ(z)`.
"""
struct MLQuad{T<:AbstractFloat}
    α::T
    β::T
    qs::QSum{T}
    sep::T
end

"""
`MLQuadPos` holds data for computing the Mittag--Leffler function `Eαβ(x)`
when `x > 0`.
"""
struct MLQuadPos{T<:AbstractFloat}
    α::T
    β::T
    qs::QSum{T}
    sep::T
end

"""
`MLQuadNeg` holds data for computing the Mittag--Leffler function `Eαβ(x)`
when `x < 0`.
"""
struct MLQuadNeg{T<:AbstractFloat}
    α::T
    β::T
    qs::QSum{T}
    sep::T
end

"""
    QSumP(T, N)

Construct quadrature sum using an optimised parabolic contour.
"""
function QSumP(::Type{T}, N) where T <: AbstractFloat
    fN = convert(T, N)
    h = 3 / fN
    μ = π * fN / 12
    w = OffsetArray{Complex{T}}(undef, 0:N)
    C = OffsetArray{Complex{T}}(undef, 0:N)
    w[0] = Complex(μ, 0)
    C[0] = Complex(exp(μ), 0)
    for n = 1:N
        one_plus_iun = Complex(1, n*h)
        w[n] = μ * one_plus_iun^2
        C[n] = exp(w[n]) * one_plus_iun
    end
    A = one(T) / 4
    return QSum(w, C, A)
end

"""
    QSumH(T, N)

Construct quadrature sum using an optimised hyperbolic contour.
"""
function QSumH(::Type{T}, N) where T <: AbstractFloat
    ϕ = parse(T, ϕ_opt)
    a, b = coef_h(ϕ), coef_μ(ϕ)
    h = a / N
    μ = b * N
    w = OffsetArray{Complex{T}}(undef, 0:N)
    C = OffsetArray{Complex{T}}(undef, 0:N)
    w0 = μ * ( 1 - sin(ϕ) )
    w[0] = Complex(w0, 0)
    for n = 0:N
        un = n * h
        iunmϕ = Complex(-ϕ, un)
        w[n] = μ * ( 1 + sin(iunmϕ) )
        C[n] = exp(w[n]) * cos(iunmϕ)
    end
    A = ( 4ϕ - π ) / 2
    return QSum(w, C, A)
end

"""
    MLQuad(α, β, N, contour, sep=0.2)

Creates an `MLQuad` object for computing `Eαβ(z)`.  The `contour` must be
either a `:parabola` or `:hyperbola`.
"""
function MLQuad(α::T, β::T, N::Integer,
                contour::Symbol, sep=one(T)/5) where T <: AbstractFloat
    if α < 0
        throw(DomainError(α, "α must be non-negative"))
    end
    if contour == :hyperbola
        qs = QSumH(T, N)
    elseif contour == :parabola
        qs = QSumP(T, N)
    else
        throw(DomainError(contour, "unrecognised"))
    end
    return MLQuad(α, β, qs, sep)
end

function (E::MLQuad{T})(z::Complex{T}) where T <: AbstractFloat
    return mlfunc(E.α, E.β, z, E.qs, E.sep)
end

function mlfunc(α::T, β::T, z::Complex{T}, 
                qs::QSum{T}, sep::T) where T <: AbstractFloat
    w, C, A = qs.w, qs.C, qs.A
    N = axes(w, 1)[end]
    m = ceil(Int64, α)
    fm = ceil(α)
    r, θ = abs(z), angle(z)
    if m == 1
        if abs(θ) > α*π
            s = C[0] * w[0]^(α-β) / ( w[0]^α - z )
            for n = 1:N
                f_wn = w[n]^(α-β) / ( w[n]^α - z ) 
                f_wnbar = conj(w[n])^(α-β) / ( conj(w[n])^α - z )
                s += C[n] * f_wn + conj(C[n]) * f_wnbar
            end
            Eαβ = A * s
        else
            γ = r^(1/α) * exp(Complex(zero(T), θ/α))
            s = C[0] * f(α, β, w[0], z, sep)
            for n = 1:N
                s += (        C[n] * f1(α, β, w[n], z, sep)
                      + conj(C[n]) * f1(α, β, conj(w[n]), z, sep) )
            end
            Eαβ = γ^(1-β) * exp(γ) / α + A * s
        end
    else
        s = zero(T)
        for k = 0:m-1
            mth_root = z^(one(T)/m)
            phase = exp(Complex(zero(T), 2π*k/fm))
            s += mlfunc(α/m, β, mth_root*phase, qs, sep)
        end
        Eαβ = s / m
    end
    return Eαβ
end

function f1(α::T, β::T, w::Complex{T}, z::Complex{T},
           sep::T) where T <: AbstractFloat
    r, θ = abs(z), angle(z)
    γ = r^(1/α) * exp(Complex(zero(T), θ/α))
    ϵ = ( w - γ ) / γ
    if abs(ϵ) > sep
        return w^(α-β) / (w^α-z) - 1 / (α*ϵ*γ^β)
    else
        return ( ψ1(α-β, ϵ) - ψ2(α, ϵ)/α ) / ( γ^β * ψ1(α, ϵ) )
    end
end

function (E::MLQuad{T})(x::T) where T <: AbstractFloat
    α, β, qs, sep = E.α, E.β, E.qs, E.sep
    if α < 2
        if x > 0
            return mlfunc_pos(α, β, x, qs, sep)
        elseif x < 0
            if α ≤ 1
                return mlfunc_neg1(α, β, -x, qs, sep)
            else
                return mlfunc_neg2(α, β, -x, qs, sep)
            end
        else
            return rΓ(β)
        end
    else
        return real(E(Complex(x, zero(T))))
    end
end

function mlfunc_pos(α::T, β::T, x::T, 
                qs::QSum{T}, sep::T) where T <: AbstractFloat
    # Compute Eαβ(x) for x≥0 and 0 < α < 2.
    w, C, A = qs.w, qs.C, qs.A
    N = axes(w, 1)[end]
    s = zero(T)
    for n = 1:N
        s += real( C[n] * f_pos(α, β, w[n], x, sep) )
    end
    s = real( C[0] * f_pos(α, β, w[0], x, sep) ) + 2s
    return ( x^((1-β)/α) * exp(x^(1/α)) / α ) + A * s
end

function f_pos(α::T, β::T, w::Complex{T}, x::T, sep::T) where T <: AbstractFloat
    xa = x^(1/α)
    ϵ = ( w - xa ) / xa
    if abs(ϵ) > sep
        return w^(α-β) / (w^α-x) - 1 / (α*ϵ*x^(β/α))
    else
        return ( ψ1(α-β, ϵ) - ψ2(α, ϵ)/α ) / ( x^(β/α) * ψ1(α, ϵ) )
    end
end

function mlfunc_neg1(α::T, β::T, x::T, 
                qs::QSum{T}, sep::T) where T <: AbstractFloat
    # Compute Eαβ(-x) for x≥0 and 0 < α < 1.
    w, C, A = qs.w, qs.C, qs.A
    N = axes(w, 1)[end]
    s = zero(T)
    for n = 1:N
        s += real( C[n] * w[n]^(α-β) / ( x + w[n]^α ) )
    end
    return A * ( real( C[0] * w[0]^(α-β) / ( x + w[0]^α ) ) + 2s )
end

function mlfunc_neg2(α::T, β::T, x::T, 
                qs::QSum{T}, sep::T) where T <: AbstractFloat
    # Compute Eαβ(-x) for x≥0 and 1 < α < 2.
    w, C, A = qs.w, qs.C, qs.A
    N = axes(w, 1)[end]
    s = zero(T)
    for n = 1:N
        f₊, f₋ = f_neg(α, β, w[n], x, sep)
        s += real( C[n] * ( f₊ + f₋ ) )
    end
    f₊, f₋ = f_neg(α, β, w[0], x, sep)
    s = real( C[0] * ( f₊ + f₋ ) ) + 2s
    sum_residues = (2/α) * ( x^((1-β)/α) * exp(x^(1/α)*cos(π/α))
             * cos( π*(1-β)/α + x^(1/α)*sinpi(1/α) ) )
    return sum_residues + A * s
end

function f_neg(α::T, β::T, w::Complex{T}, x::T, sep::T) where T <: AbstractFloat
    γ₊ = x^(1/α) * exp(complex(zero(T), π/α))
    ϵ₊ = ( w - γ₊ ) / γ₊
    γ₋= x^(1/α) * exp(complex(zero(T), -π/α))
    ϵ₋ = ( w - γ₋) / γ₋
    if abs(ϵ₊) < sep
        numer  = ( ( w - γ₋ ) * ( ψ1(α-β, ϵ₊) - ψ2(α, ϵ₊)/α )
                  - γ₊ * ψ1(α, ϵ₊)/α ) 
        denom = γ₊^β * ψ1(α, ϵ₊) * ( w - γ₋ + γ₊*ϵ₊ )  
        f₊ = numer / denom
        f₋ = ( γ₊^(1-β) * (1+ϵ₊)^(α-β) / ( ψ1(α,ϵ₊)*( 2w - γ₊ - γ₋ ) )
              - γ₋^(1-β) / ( α * ( w - γ₋ ) ) )
    elseif abs(ϵ₋) < sep
        numer  = ( ( w - γ₊ ) * ( ψ1(α-β, ϵ₋) - ψ2(α, ϵ₋)/α ) 
                  - ψ1(α, ϵ₋) * γ₋ / α )
        denom = γ₋^β * ψ1(α, ϵ₋) * ( 2w - γ₋ - γ₊ )  
        f₋ = numer / denom
        f₊ = ( γ₋^(1-β) * (1+ϵ₋)^(α-β) / ( ψ1(α,ϵ₋)*( 2w - γ₊ - γ₋ ) )
              - γ₊^(1-β) / ( α * ( w - γ₊ ) ) )
    else
        f₊ = ( w^(α-β) * ( w - γ₋ ) / ( ( w^α + x ) * ( 2w - γ₊ - γ₋ ) )
              - γ₊^(1-β) / ( α * ( w - γ₊ ) ) )
        f₋ = ( w^(α-β) * ( w - γ₊ ) / ( ( w^α + x ) * ( 2w - γ₋ - γ₊ ) )
              - γ₋^(1-β) / ( α * ( w - γ₋ ) ) )
    end
    return f₊, f₋
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

end # module
