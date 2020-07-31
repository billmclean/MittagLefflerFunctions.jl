import MittagLefflerFunctions.Quadrature.QSumP
import MittagLefflerFunctions.Quadrature.QSumH
using PyPlot

figure(1)
marker = ("^", "o", "s", "v")
N = (5, 10, 15, 20)
for k = 1:4
    qs = QSumP(Float64, N[k])
    plot(real.(qs.w), imag.(qs.w), marker[k], label="N = $(N[k])")
end
legend()
grid(true)
axis((-50, 10, -5, 50))
xlabel(L"$\Re w$", fontsize=12)
ylabel(L"$\Im w$", fontsize=12)
savefig("figure4L.pdf")

figure(2)
for k = 1:4
    qs = QSumH(Float64, N[k])
    plot(real.(qs.w), imag.(qs.w), marker[k], label="N = $(N[k])")
end
legend()
grid(true)
axis((-50, 10,-5, 50))
xlabel(L"$\Re w$", fontsize=12)
ylabel(L"$\Im w$", fontsize=12)
savefig("figure4R.pdf")

figure(4)
for k = 1:4
    qs = QSumH(Float64, N[k])
    semilogy(collect(0:N[k]), abs.(qs.C), marker[k], label="N = $N")
end
legend()
grid(true)
xlabel(L"$n$", fontsize=12)
ylabel(L"$|C_n|$", fontsize=12)
xylims = axis()
savefig("figure5L.pdf")

figure(3)
for k = 1:4
    qs = QSumP(Float64, N[k])
    semilogy(collect(0:N[k]), abs.(qs.C), marker[k], label="N = $N")
end
legend()
grid(true)
#axis((-1.0, 21.0, 7.74879598668249e-22, 6034.9888883983585))
axis(xylims)
xlabel(L"$n$", fontsize=12)
ylabel(L"$|C_n|$", fontsize=12)
savefig("figure5R.pdf")
