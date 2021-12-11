using OffsetArrays
using SpecialFunctions: gamma, erfcx
using PyPlot
using Polynomials
using Printf
using LinearAlgebra

Γ(z) = gamma(z)

function pade_matrix(α::T, β::T, 
                     m::Integer, n::Integer) where T <: AbstractFloat
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
    A = zeros(T, 2r, 2r+1)
    if m ≥ r+1
        row = 0
        for k = 0:r-1
            row += 1
            A[row,k+1] = one(T)
            for j = 0:k
                A[row,r+1+j] = -a[k-j]
            end
        end
        for k = r:m-1
            row += 1
            for j = 0:r
                A[row,r+1+j] = -a[k-j]
            end
        end
        for k = r-n+1:r-1
            row += 1
            A[row,k+1] = one(T)
            for j = k+1:r
                A[row,r+1+j] = -b[j-k]
            end
        end
    else
        row = 0
        for k = 0:m-1
            row += 1
            A[row,k+1] = one(T)
            for j = 0:k
                A[row,r+1+j] = -a[k-j]
            end
        end
        for k = -(r-m):-1
            row += 1 
            for j = 0:r
                A[row,r+1+j] = -b[j-k]
            end
        end
        for k = 0:r-1
            row += 1
            A[row,k+1] = one(T)
            for j = k+1:r
                A[row,r+1+j] = -b[j-k]
            end
        end
    end
    return A
end

function pade_svd(α::T, β::T, m::Integer, n::Integer) where T <: AbstractFloat
    A = pade_matrix(α, β, m, n)
    p, q = pade_svd(A)
    return p, q
end

function pade_svd(A::Matrix{T}) where T <: AbstractFloat
    r = div(size(A,1), 2)
    F = svd(A, full=true)
    p = F.V[1:r,2r+1]
    q = F.V[r+1:2r+1,2r+1]
    q1 = q[1]
    p ./= q1
    q ./= q1
    return p, q
end

function pade_lu(A::Matrix{T}) where T <: AbstractFloat
    r = div(size(A,1), 2)
    F = lu(A)
    y = zeros(T, 2r+1)
    y[2r+1] = one(T)
    for i = 2r:-1:1
        s = zero(T)
        for j = i+1:2r+1
            s += F.U[i,j] * y[j]
        end
        y[i] = -s / F.U[i,i]
    end
    p = y[1:r]
    q = y[r+1:2r+1]
    μ = q[1]
    for k = 1:r
        p[k] /= μ
    end
    for k = 1:r+1
        q[k] /= μ
    end
    return p, q, y
end

#T = BigFloat
T = Float64
α = 1 / parse(T, "2")
β = parse(T, "1")

x = range(zero(T), parse(T, "10"), length=201)
y = erfcx.(x)
figure(1)
r = 3
for m = r:5
    n = 2r + 1 - m
    A = pade_matrix(α, β, m, n)
    p, q = pade_lu(A)
    P = Polynomial(p)
    Q = Polynomial(q)
    err = P.(x) ./ Q.(x) - y
    plot(x, err, label="m, n = $m, $n")
end
legend()
grid(true)

figure(2)
for (r, m, n) in [(4, 5, 4), (8, 9, 8), (12, 13, 12)]
    A = pade_matrix(α, β, m, n)
    p, q = pade_lu(A)
    P = Polynomial(p)
    Q = Polynomial(q)
    err = P.(x) ./ Q.(x) - y
    semilogy(x, abs.(err), label="m, n = $m, $n")
end
legend()
grid(true)

figure(3)
p, q = pade_svd(1.85, 1.0, 12, 11)
P = Polynomial(p)
Q = Polynomial(q)
R(x) = P(x) / Q(x)
plot(x, R.(x))
grid(true)
