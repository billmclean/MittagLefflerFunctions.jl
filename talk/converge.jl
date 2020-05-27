using ArgCheck
using Printf

A(ϕ) = acosh(2ϕ/((4ϕ-π)*sin(ϕ)))

function Q(x::T, α::T, N::Integer) where T <: AbstractFloat
    @argcheck 0 < α < 1
    ϕ = parse(T, "1.17210")
    h = A(ϕ) / N
    μ = ( (4π*ϕ-π^2) / A(ϕ) ) * N
    w0 = μ * ( 1 - sin(ϕ) )
    c0 = exp(w0) * w0^(α-1) * μ * cos(ϕ) 
    a0 = w0^α
    s = (c0/2) / ( x + a0 )
    for n = 1:N
        un = n * h
        wn = Complex(μ * ( 1 - cosh(un) * sin(ϕ) ), μ * sinh(un) * cos(ϕ))
        wnα = wn^α
        an = real(wnα)
        bn = imag(wnα)
        term = exp(wn) * (wnα/wn) * μ * cos(Complex(-ϕ, un)) 
        cn = real(term)
        dn = imag(term)
        s += ( cn*(x+an) + bn*dn) / ( (x+an)^2 +bn^2 )
    end
    return h * s / π
end

T = BigFloat
ten = parse(T, "10.0")
α = 3/ten
for N = 4:4:20
    @printf("%4d & ", N)
    err = Q(zero(T), α, N) - 1
    @printf("%10.2e \\\\\n", err)
end


