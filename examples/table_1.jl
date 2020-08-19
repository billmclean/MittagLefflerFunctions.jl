import MittagLefflerFunctions: rΓ, MLQuad
using Printf
import SpecialFunctions: gamma

Γ(x) = gamma(x)

α = 0.7
β = 1.0
tol = 1e-12
Eref = MLQuad(α, β, 14, :hyperbola)
xvals = collect(5.0:10.0:55.0)

function asymp_coeffs(α, β, N)
    s = 0.0
    σ = zeros(N)
    τ = zeros(N)
    sgn = 1.0
    s = 0.0
    for n = 1:N
        if n*α < β
            σ[n] = 1.0
            τ[n] = 1 / Γ(β-n*α)
        else
            σ[n] = -sin(π*(n*α-β))
            τ[n] = Γ(1+n*α-β) / π
        end
    end
    return σ, τ
end

function table_row(α, β, x, tol)
    Nopt = x^(1/α) / α
    N = ceil(Int64, Nopt)
    σ, τ = asymp_coeffs(α, β, N)
    z = -x
    powz = 1.0
    m = 0
    next = NaN
    last = NaN
    s = 0.0
    for n = 1:N-1
        powz *= z
        s -= σ[n] * τ[n] / powz
        last = τ[n] / abs(powz)
        m = n + 1
        if last < tol
            break
        end
    end
    powz *= z
    next = τ[m] / abs(powz)
    return s, m, Nopt, last, next
end

@printf("α = %0.2f, β = %0.2f, tol = %0.2g\n\n", α, β, tol)
@printf("%6s  %5s  %8s  %10s  %10s  %10s\n\n", 
        "x", "m", "Nx", "error", "m-1", "m")
for x in xvals
    s, m, Nopt, last, next = table_row(α, β, x, tol)
    err = Eref(-x) - s
    @printf("%6.2g& %5d& %8.1f& %10.2e& %10.2e& %10.2e\\\\\n", 
            x, m, Nopt, err, last, next)
end
