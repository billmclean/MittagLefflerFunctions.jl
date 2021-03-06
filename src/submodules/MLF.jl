"""
A direct translation of the Matlab function mlf by Igor Podlubny and Martin 
Kacenak.  https://www.mathworks.com/matlabcentral/fileexchange/8738-mittag-leffler-function
"""
module MLF

using ..MittagLefflerFunctions: rΓ

export mlf

"""
    mlf(α, β, z, p)

Evaluate the Mittag-Leffler function E_αβ(z) with accuracy 10^(-p).
"""
function mlf(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    ten = parse(T, "10.0")
    r0 = radius(α, β, z, p)
    if α == one(T) && β == one(T)
        E = exp(z)
    else
        if (   ( α < 1 && abs(z) ≤ 1 ) 
            || ( ( 1 ≤ α < 2 ) && abs(z) ≤ floor(20/(2.1-α)^(5.5-2α))) 
            || (  α ≥ 2 && abs(z) ≤ 50 ) )
            E = mlf1(α, β, z, p)
        else
            if  α<=1 && abs(z) ≤ floor(5*α+ten) 
                if (   ( abs(angle(z)) > π*α ) 
                    && ( abs(abs(angle(z)) - (π*α)) > ten^(-p) ) )
                    if β ≤ 1
                        E = mlf2(α, β, z, p, r0)
                    else
                        E = mlf3(α, β, z, p, r0)
                    end
                elseif (   abs(angle(z)) < π*α 
                        && abs(abs(angle(z))-(π*α)) > ten^(-p) )
                    if β <= 1
                        E = mlf4(α, β, z, p, r0)
                    else
                        E = mlf5(α, β, z, p, r0)
                    end
                else
                    E = mlf6(α, β, z, p, r0)
                end
            else
                if α ≤ 1
                    if abs(angle(z)) < (π*α/2+min(π,π*α))/2
                        E = mlf7(α, β, z, p)
                    else
                        E = mlf8(α, β, z, p)
                    end
                else
                    if α ≥ 2
                        E = mlf9(α, β, z, p)
                    else
                        E = mlf10(α, β, z, p)
                    end
                end
            end
        end
    end
    return E
end

function mlf(α::T, β::T, x::T, p::Integer) where T <: AbstractFloat
    z = Complex(x, zero(T))
    return real(mlf(α, β, z, p))
end

function radius(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    ten = parse(T, "10.0")
    if β < 0 
        rc = ( -2*log( ten^(-p)*π / (6*(abs(β)+2)*(2*abs(β))^(abs(β))) ) )^α
    else  
        rc = ( -2*log( ten^(-p)*π / 6 ) )^α
    end
    r0 = max(one(T), 2*abs(z), rc)
    return r0
end

function mlf1(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    oldsum=0
    k = 0
    powz = Complex(one(T),zero(T))
    while β + k*α ≤ 0 
        powz *= z # = z^k
        k += 1
    end
    newsum = powz * rΓ(β+k*α)
    while newsum != oldsum
        oldsum = newsum
        k += 1
        powz *= z
        term = powz * rΓ(β+k*α)
        newsum += term
        k += 1
        powz *= z
        term = powz * rΓ(β+k*α)
        newsum += term
    end
    E = newsum
    return E
end

function mlf2(α::T, β::T, z::Complex{T}, p::Integer,
              r0::T) where T <: AbstractFloat
    E = Romberg(zero(T), r0, p) do r
        return K(r, α, β, z)
    end
    return E
end

function mlf3(α::T, β::T, z::Complex{T}, p::Integer,
              r0::T) where T <: AbstractFloat
    ϵ = one(T)
    E1 = Romberg(ϵ, r0, p) do r
        return K(r, α, β, z) 
    end
    E2 = Romberg(-π*α, π*α, p) do r
        return P(r, α, β, z, ϵ) 
    end
    E = E1 + E2
    return E
end

function mlf4(α::T, β::T, z::Complex{T}, p::Integer,
              r0::T) where T <: AbstractFloat
    E1 = Romberg(zero(T), r0, p) do r
        return K(r, α, β, z)
    end
    E = E1 + z^((1-β)/α) * exp(z^(1/α)) / α 
    return E
end

function mlf5(α::T, β::T, z::Complex{T}, p::Integer,
              r0::T) where T <: AbstractFloat
    ϵ = abs(z)/2
    E1 = Romberg(ϵ, r0, p) do r
        return K(r, α, β, z)
    end
    E2 = Romberg(-π*α, π*α, p) do r
        return P(r, α, β, z, ϵ)
    end
    E = E1 + E2 + (z^((1-β)/α)) * (exp(z^(1/α))/α)
    return E
end

function mlf6(α::T, β::T, z::Complex{T}, p::Integer,
              r0::T) where T <: AbstractFloat
    ϵ = abs(z) + one(T)/2
    E1 = Romberg(ϵ, r0, p) do r
        return K(r, α, β, z)
    end
    E2 = Romberg(-π*α, π*α, p) do r
        return P(r, α, β, z, ϵ) 
    end
    E = E1 + E2
    return E
end

function mlf7(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    newsum = z^((1-β)/α) * exp(z^(1/α))/α
    powz = Complex(one(T), zero(T))
    for k=1:floor(Int64, p/log10(abs(z)))
        powz /= z # = z^(-k)
        newsum -= powz * rΓ(β-k*α)
    end
    E = newsum
    return E
end

function mlf8(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    newsum = zero(Complex{T})
    powz = Complex(one(T), zero(T))
    for k=1:floor(Int64, p/log10(abs(z)))
        powz /= z # = z^(-k)
        newsum -= powz * rΓ(β-k*α)
    end
    E = newsum
    return E
end

function mlf9(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    m = floor(Int64, α/2)
    s = zero(Complex{T})
    for h = 0:m
        zn = z^(one(T)/(m+1)) * exp_i(2π*h/(m+1))
        s += mlf(α/(m+1), β, zn, p)
    end
    E = (one(T)/(m+1)) * s
    return E
end

function mlf10(α::T, β::T, z::Complex{T}, p::Integer) where T <: AbstractFloat
    E = ( mlf(α/2, β, sqrt(z), p)
        + mlf(α/2, β, -sqrt(z), p) ) / 2
    return E
end

function exp_i(t::T) where T <: AbstractFloat
    return exp(Complex(zero(T),t))
end

function K(r::T, α::T, β::T, z::Complex{T}) where T <: AbstractFloat
    return r^((1-β)/α) * exp(-r^(1/α)) * ( (r*sinpi(1-β)-z*sinpi(1-β+α))
                                          /(π*α*(r^2-2*r*z*cospi(α)+z^2)) )
end

function P(r::T, α::T, β::T, z::Complex{T}, ϵ::T) where T <: AbstractFloat
    w = ϵ^(1/α) * sin(r/α) + r * ( 1 + (1-β)/α )
    numer =  exp((ϵ^(1/α)) * cos(r/α)) * exp_i(w)
    return ( ϵ^(1+(1-β)/α) / (2π*α) ) * numer / ( ϵ*exp_i(r) - z )
end

function Romberg(f::Function, a::T, b::T, order=6) where T <: AbstractFloat
    R = zeros(Complex{T}, 2, order)
    h = b - a
    R[1,1] = ( h/2 ) * ( f(a) + f(b) )
    ipower = 1
    for i = 2:order
        s = zero(Complex{T})
        for j = 1:ipower
            s += f( a + (2j-1)*h/2 )
        end
        R[2,1] = ( R[1,1] + h * s ) / 2
        pow4 = one(T)
        for k = 1:i-1
            pow4 *= 4
            R[2,k+1] = ( pow4 * R[2,k] - R[1,k] ) / ( pow4 - 1 )
        end
        for j = 0:i-1
            R[1,j+1] = R[2,j+1]
        end
        ipower *= 2
        h /= 2
    end
    return R[1,order]
end

end # module
