module Quadrature

using OffsetArrays

export MLQuad1

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

end # module
