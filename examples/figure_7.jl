using MittagLefflerFunctions
using PyPlot
using SpecialFunctions: erfcx

α = 1/2
β = 1.0

#x = range(0, 10, length=201)
x = range(-5, 3, length=300)
y = erfcx.(-x)

Nvals = collect(4:2:10)

function draw_errors(contour::Symbol)
    for row = 1:length(Nvals)
        N = Nvals[row]
        E = MLQuad(α, β, N, contour)
        err = abs.( y - E.(x) )
        semilogy(x, err, label="N = $N") 
    end
    legend(loc="upper right")
    grid(true)
    xlabel(L"$x$")
end

figure(1)
draw_errors(:parabola)
xylims = (-5.4, 3.4, 3e-12, 6e-4)
axis(xylims)
savefig("figure7L.pdf")

figure(2)
draw_errors(:hyperbola)
savefig("figure7R.pdf")
axis(xylims)

