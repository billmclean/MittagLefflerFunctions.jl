using PyPlot
using OffsetArrays
using SpecialFunctions: erfcx, gamma

const ϕ = 1.172104
const a = 1.081792
const b = 4.492075

α = 0.5
β = 1.0
F(w, x, α, β) = w^(α-β) / (w^α+x)
F0(w, x, α, β) = -x / (w^β*(w^α+x))

function contour(u, N)
    μ = b * N
    iumϕ = Complex(-ϕ, u)
    w = μ * ( 1 + sin(iumϕ) )
    dw = μ * cos(iumϕ)
    return w, dw
end

function points(N)
    w = OffsetArray{Complex{Float64}}(undef, 0:N)
    dw = OffsetArray{Complex{Float64}}(undef, 0:N)
    h = a / N
    for n = 0:N
        w[n], dw[n] = contour(n*h, N)
    end
    return w, dw, h
end

function Taylor(x, α, β, N)
    s = 1.0 / gamma(β)
    powx = 1.0 
    for n = 1:N
        powx *= -x
        s += powx / gamma(β+n*α)
    end
    return s
end

function newI(x, α, β, N)
    w, dw, h = points(N)
    s = exp(w[0]) * F(w[0], x, α, β) * dw[0]
    for n = 1:N
        s += 2 * exp(w[n]) * F(w[n], x, α, β) * dw[n]
    end
    return h * real(s) / (2π)
end

figure(1)
for N in [5, 10, 15, 20]
    w, dw, h = points(N)
    plot(real.(w), imag.(w), "o", markersize=2, label="N = $N")
end
legend()
xlabel(L"$\Re w$")
ylabel(L"$\Im w$")
grid(true)
savefig("points.pdf")

figure(2)
x = 2.0
N = 15
w, dw, h = points(N)
u = Float64[ n*h for n = 0:N ]
integrand = F.(w, x, α, β) .* exp.(w) .* dw / Complex(0.0, 2π)
plot(u, real.(integrand[0:N]), "+-", label="real part")
#plot(u, imag.(integrand[0:N]), "x-", label="imag part")
#legend()
grid(true)
xlabel(L"$u$")
title("N = $N")
savefig("integrand.pdf")
#    semilogy(u, abs.(real.(integrand[0:N])),
#             u, abs.(imag.(integrand[0:N])))

figure(3)
x = range(0.0, 5.0, length=201)
nn = [ n for n=1:length(x) if x[n] ≤ 0.2 ]
y1 = erfcx.(x)
N = 10
y2 = newI.(x, 0.5, 1.0, N)
plot(x, y1-y2)
#y3 = Taylor.(x[nn], 0.5, 1.0, N)
#plot(x, y1-y2, x[nn], y1[nn]-y3)
#plot(x, y2, x[nn], y3)
xlabel(L"$x$")
title("N = $N")
grid(true)
savefig("error.pdf")
