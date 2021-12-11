using OffsetArrays
using SpecialFunctions: gamma, erfcx
using Polynomials
using PyPlot
using Printf

T = Float64
α = 1 / parse(T, "2")
β = one(T)
a = OffsetArray{T}(undef, 0:4)
b = OffsetArray{T}(undef, 0:4)
a[0] = 1 / gamma(β)
b[0] = zero(T)
let
pow = one(T)
for k = 1:4
    pow = -pow
    a[k] = pow / gamma(β+k*α)
    b[k] = pow * sinpi(k*α-β) * gamma(1+k*α-β) / π
end
end

A_43 = [ 1 0 -a[0]    0     0
      0 1 -a[1] -a[0]    0
      0 0 -a[2] -a[1] -a[0]
      1 0 -b[0] -b[1] -b[2]
      0 1    0  -b[0] -b[1] ]

rhs_43 = [ a[1], a[2], a[3], 0, 0 ]
pq = A_43 \ rhs_43

p = OffsetArray{T}(undef,0:2)
q = OffsetArray{T}(undef,0:3)
p[0] = a[0]
p[1:2] = pq[1:2]
q[0] = one(T)
q[1:3] = pq[3:5]

P_43 = Polynomial(p[0:2])
Q_43 = Polynomial(q[0:3])
R_43(x) = P_43(x) / Q_43(x)

A_52 = [ 1 0 -a[0]    0     0
         0 1 -a[1] -a[0]    0
         0 0 -a[2] -a[1] -a[0]
         0 0 -a[3] -a[2] -a[1]
         0 1    0  -b[0] -b[1] ]

rhs_52 = [ a[1], a[2], a[3], a[4], 0 ]
pq = A_52 \ rhs_52

p = OffsetArray{T}(undef,0:3)
q = OffsetArray{T}(undef,0:3)
p[0] = a[0]
p[1:2] = pq[1:2]
p[3] = zero(T)
q[0] = one(T)
q[1:3] = pq[3:5]

P_52 = Polynomial(p[0:2])
Q_52 = Polynomial(q[0:3])
R_52(x) = P_52(x) / Q_52(x)

figure(1)
x = range(0, 10, length=201)
Eα = erfcx.(x)
plot(x, R_43.(x)-Eα, x, R_52.(x)-Eα)
grid(true)
legend((L"$m=4, n=3$", L"m=5, n=2"))
xlabel(L"$x$", fontsize=12)
savefig("pade_43_52.pdf")

for k = 1:3
    @printf("%4d& %8.5f& %8.5f& %8.5f& %8.5f\\\\\n", k-1,
            P_43.coeffs[k], Q_43.coeffs[k], P_52.coeffs[k], Q_52.coeffs[k])
end
@printf("%4d& %8.5f& %8.5f& %8.5f& %8.5f\n", 3,
        0.0, Q_43.coeffs[4], 0.0, Q_52.coeffs[4])
