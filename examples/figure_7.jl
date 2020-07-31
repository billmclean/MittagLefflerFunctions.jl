using MittagLefflerFunctions
using PyPlot
using SpecialFunctions: erfcx

α = 1/2
β = 1.0

#x = range(0, 10, length=201)
x = range(-5, 3, length=201)
y = erfcx.(-x)

Nvals = collect(4:2:10)

function draw_errors(contour::Symbol)
    for row = 1:length(Nvals)
        N = Nvals[row]
        E = MLQuad(α, β, N, contour)
        err = abs.( y - E.(x) )
        semilogy(x[2:end], err[2:end], label="N = $N") # omit x=0
    end
legend(loc="upper right")
grid(true)
xlabel(L"$x$")
end

figure(2)
draw_errors(:hyperbola)
savefig("figure7R.pdf")
xylims = axis()

figure(1)
draw_errors(:parabola)
axis(xylims)
savefig("figure7L.pdf")
