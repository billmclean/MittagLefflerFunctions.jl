module MittagLefflerFunctions

import SpecialFunctions
import GaussQuadrature
using OffsetArrays
using ArgCheck

export chebyshev_coefs!, chebyshev_sum

Γ = SpecialFunctions.gamma

"""
    ML_power_series(α, β, x, nmax, tol)

Sum the first `nmax` terms of the Taylor series expansion of `E_αβ(x)`,
or stop when the next term is smaller than `tol`.
"""
function ML_power_series(α::T, β::T, x::T, nmax::Integer,
                            tol::T) where T <: AbstractFloat
    if α ≤ 0 || β ≤ 0
        throw(ArgumentError("require α > 0  and β > 0"))
    end
    Eαβ = one(T) / Γ(β)
    powx = one(T)
    for n = 1:nmax
        powx *= x
        term = powx / Γ(β+n*α)
        Eαβ += term
        if abs(term) ≤ tol
            break
        end
    end
    return Eαβ
end

"""
    MLneg_asymptotic_series(α, β, x, tol, nmax)

Sum the first `nmax` terms of the asymptotic expansion of `E_α(-x)`,
or stop when the next term is smaller than `tol`.
"""
function MLneg_asymptotic_series(α::T, β::T, x::T, nmax::Integer
                                ) where T <: AbstractFloat
    if α ≤ 0 || α ≥ 2 || β ≤ 0
        throw(ArgumentError("require 0 < α < 2 and β > 0"))
    end
    recipx = one(T) / x
    powrx = recipx
    βmα = β - α
    if βmα > 0
        c1 = 1 / Γ(βmα)
    else
        c1 = Γ(1-βmα) * sinpi(βmα) / π
    end
    Eαβ = c1 * powrx 
    for n = 2:nmax
        βmnα = β - n * α
        if βmnα > 0
            cn = 1 / Γ(βmn*α)
        else
            cn = Γ(1-βmnα) * sinpi(βmnα) / π
        end
        powrx *= -recipx
        Eαβ += cn * powrx
    end
    return Eαβ
end

"""
    MLpos_asymptotic_series(α, β, x, tol, nmax)

Sum the first `nmax` terms of the asymptotic expansion of `E_α(+x)`,
or stop when the next term is smaller than `tol`.
"""
function MLpos_asymptotic_series(α::T, β::T, x::T, nmax::Integer
                            ) where T <: AbstractFloat
    if α ≤ 0 || α ≥ 2 || β ≤ 0
        throw(ArgumentError("require 0 < α < 2 and β > 0"))
    end
    Eαβ = ( x^((1-β)/α) * exp(x^(1/α)) / α 
          + MLneg_asymptotic_series(α, β, -x, nmax) )
    return Eαβ
end

"""
    MLpos_asymptotic_C(α, β, N)

Constant used to estimate the remainder term in the asymptotic expansion
of `Eαβ(x)`.
"""
function MLpos_asymptotic_C(α::T, β::T, N::Integer) where T <: AbstractFloat
    a = sinpi( β - (N+2)*α )
    b = sinpi(α)
    c = cospi(α)
    if abs(a) < eps(A)
        if c ≥ 0
            Cplus = 1 / abs(b)
        else
            Cplus = abs(b)
        end
    else
        # solve ay^2 + 2by - (a+2bc) = 0
        A = b / a
        B = 1 + 2A * c
        # y^2 + 2Ay - B = 0 so (y+A)^2 = B + A^2
        if B + A^2 < 0 # no real solutions
            Cplus = abs(b)
        else
            D = sqrt(B+A^2) # so y = -A ± D
            if A ≥ 0
                yminus = -A - D
                yplus = -B / yminus
            else
                yplus = -A + D
                yminus = -B / yplus
            end
            Cplus = abs(b)
            if yplus > 0
                Cplus = max(Cplus, abs((a*yplus+b)/((yplus-c)^2+b^2)))
            end
            if yminus > 0
                Cplus = max(Cplus, abs((a*yminus+b)/((yminus-c)^2+b^2)))
            end
        end
    end
    return Cplus
end

function MLneg_integral(α::T, β::T, λ::T, t::T, nmax::Integer
                        ) where T <: AbstractFloat
    F(z) = 1 / ( z^β + λ * z^(β-α) )
    ϕ = parse(T, "1.1721")
    h = parse(T, "1.0818") / nmax
    μ = parse(T, "4.4921") * nmax / t
    z0 = μ * ( 1 - sin(ϕ) )
    s = exp(z0*t) * F(z0) * cos(ϕ) / 2
    for n = 1:nmax
        iunmϕ = Complex(-ϕ, n*h)
        zn = μ * ( 1 + sin(iunmϕ) )
        s += real( exp(zn*t) * F(zn) * cos(iunmϕ) )
    end
    return μ * t^(1-β) * h * s / π
end

function chebyshev_polys!(Cheb::OffsetArray{T}, x::T
                         ) where T <: AbstractFloat
    nmax = length(Cheb) - 1
    Cheb[0] = one(T)
    if nmax ≥ 1
        Cheb[1] = x
    end
    if nmax ≥ 2
        for n = 1:nmax-1
            Cheb[n+1] = 2x * Cheb[n] - Cheb[n-1]
        end
    end
end

function chebyshev_polys!(Cheb::OffsetArray{T}, x::AbstractVector{T}
                         ) where T <: AbstractFloat
    nmax, M = size(Cheb)
    nmax -= 1
    @argcheck Cheb.offsets == (-1, 0)
    @argcheck length(x) == M
    for m = 1:M
        Cheb[0,m] = one(T)
    end
    if nmax ≥ 1
        for m = 1:M
            Cheb[1,m] = x[m]
        end
    end
    if nmax ≥ 2
        for m = 1:M, n = 1:nmax-1
            Cheb[n+1,m] = 2x[m] * Cheb[n,m] - Cheb[n-1,m]
        end
    end
end

function chebyshev_coefs!(a::OffsetArray{T}, f::Function, M::Integer
                         ) where T <: AbstractFloat
    nmax = length(a) - 1
    @argcheck a.offsets == (-1,)
    x, w = GaussQuadrature.chebyshev(T, M)
    Cheb = OffsetArray{T}(undef, 0:nmax, 1:M)
    chebyshev_polys!(Cheb, x)
    for n = 0:nmax
        s = zero(T)
        for m = 1:M
            s += w[m] * f(x[m]) * Cheb[n,m]
        end
        a[n] = (2/π) * s
    end
end

function chebyshev_sum(a::OffsetArray{T}, x::T) where T <: AbstractFloat
    nmax = length(a) - 1
    bn = bnp1 = bnp2 = zero(T)
    for n = nmax:-1:1
        bn = a[n] + 2x*bnp1 - bnp2
        bnp2 = bnp1
        bnp1 = bn
    end
    bn = a[0] + 2x*bnp1 - bnp2
    return ( bn - bnp2 ) / 2
end

end # module
