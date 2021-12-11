using OffsetArrays
using SpecialFunctions: gamma, erfcx
using PyPlot
using Polynomials
using Printf
using LinearAlgebra

Γ(z) = gamma(z)

function pade(α::T, β::T, m::Integer, n::Integer) where T <: AbstractFloat
    if mod(m+n,2) == 0
        throw(ArgumentError("m+n must be odd"))
    end
    r = div(m+n-1,2)
    a = OffsetArray{T}(undef, 0:m-1)
    b = OffsetArray{T}(undef, 0:n-1)
    pow = one(T)
    for k = 0:m-1
        a[k] = pow / Γ(β+k*α)
        pow = -pow
    end
    b[0] = zero(T)
    pow = one(T)
    for k = 1:n-1
        pow = -pow
        b[k] = pow * sinpi(k*α-β) * Γ(1+k*α-β) / π
    end
    A = zeros(T, 2r-1, 2r-1)
    rhs = zeros(T, 2r-1)
    if m ≥ r+1
        for k = 1:r-1
            A[k,k] = one(T)
            for j = 1:k
                A[k,r-1+j] = -a[k-j]
            end
            rhs[k] = a[k]
        end
        for k = r:m-1
            for j = 1:r
                A[k,r-1+j] = -a[k-j]
            end
            rhs[k] = a[k]
        end
        for k = r-n+1:r-1
            row = k+r
            A[row,k] = one(T)
            for j = k+1:r
                col = j + r - 1
                A[row,col] = -b[j-k]
            end
        end
    else
        for k = 1:m-1
            A[k,k] = one(T)
            for j = 1:k
                col = r - 1 + j
                A[k,col] = -a[k-j]
            end
            rhs[k] = a[k]
        end
        for k = 1:n-r-1
            row = k + m - 1
            for j = 1:r
                col = r - 1 + j
                A[row,col] = -b[j+k]
            end
            rhs[row] = b[k]
        end
        for j = 1:r
            col = r - 1 + j 
            A[r,col] = -b[j]
        end
        rhs[r] = -a[0]
        for k = 1:r-1
            row = r + k
            A[row,k] = one(T)
            for j = k+1:r
                col = r - 1 + j
                A[row,col] = -b[j-k]
            end
        end
    end
    coeffs = A \ rhs
    p = OffsetArray{T}(undef, 0:r-1)
    q = OffsetArray{T}(undef, 0:r)
    p[0] = a[0]
    for k = 1:r-1
        p[k] = coeffs[k]
    end
    q[0] = one(T)
    for k = 1:r
        q[k] = coeffs[k+r-1]
    end
    return p, q, A, rhs
end

#T = BigFloat
T = Float64
r = 3
α = 1 / parse(T, "2")
β = parse(T, "1")
x = range(zero(T), parse(T, "10"), length=201)
y = erfcx.(x)
figure(1)
for m = r:5
    n = 2r + 1 - m
    @printf("%3d  %3d", m, n)
    p, q, A, rhs = pade(α, β, m, n)
    @printf("  %5.2f\n", cond(A))
    P = Polynomial(p[0:r-1])
    Q = Polynomial(q[0:r])
    err = P.(x) ./ Q.(x) - y
    plot(x, err, label="m, n = $m, $n")
end
legend()
grid(true)

figure(2)
#for (r, m, n) in [(4, 5, 4), (8, 9, 8), (12, 13, 12), (16, 17, 16)]
for (r, m, n) in [(4, 5, 4), (8, 9, 8), (12, 13, 12)]
    p, q, A, rhs = pade(α, β, m, n)
    println("cond(A) = ", cond(A))
    display(p)
    display(q)
    P = Polynomial(p[0:r-1])
    Q = Polynomial(q[0:r])
    err = P.(x) ./ Q.(x) - y
    semilogy(x, abs.(err), label="m, n = $m, $n")
end
legend()
grid(true)
