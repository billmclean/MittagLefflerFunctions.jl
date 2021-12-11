using SpecialFunctions: loggamma, gamma
using Printf: @printf

T = Float64

α = 0.6
β = 0.8
x = 20.0
tol = 1e-16
x_vals = T[ n*5 for n = 5:10 ]

function asymp_terms(x, α, β, tol)
    Nopt = ceil(Integer, x^(1/α))
    @printf("α = %0.4f, β = %0.4f\n", α, β)
    @printf("x = %0.4f, Nopt = %0d\n", x, Nopt)
    sgn = one(T)
    s = zero(T)
    @printf("\n%5s  %12s  %12s  %15s\n\n", 
            "n", "B_n/x^n", "term", "partial sum")
    for n = 1:Nopt
        A = sinpi(n*α-β)
        B = loggamma(n*α-β+1)
        Bxn = exp(B-n*log(x))
        term = sgn * A * Bxn / π
        s += term
        sgn = -sgn
        @printf("%5d  %12.4e  %12.4e  %20.15f\n", n, Bxn, term, s)
        if  Bxn < tol
            break
        end
    end
end

function minimal_x_N(α, β, tol, minx=10, maxN=100)
    first_n = 1
    huge = convert(T, typemax(1))
    for k = minx:maxN
        x = convert(T, k)
        if log(x) < floor(α * log(huge))
            Nopt = min(maxN, ceil(Integer, x^(1/α)))
        else
            Nopt = maxN
        end
        for n = first_n:Nopt 
            bound = loggamma(n*α-β+1) - n*log(x)
            @printf("%5d  %8.3f  %12.4e\n", n, x, exp(bound))
            if bound < log(tol)
                @goto done
            end
        end
        first_n = Nopt
    end
    @label done
end

asymp_terms(x, α, β, tol)

