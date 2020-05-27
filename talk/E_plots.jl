using PyPlot
using SpecialFunctions

figure(1)
x1 = range(-3, 1.8, length=201)
xhalf = range(-3, 1, length=201)
x0 = range(-3, 0.8, length=201)
plot(x1, exp.(x1), xhalf, erfcx.(-xhalf), x0, 1 ./(1 .- x0))
legend((L"$E_1$", L"$E_{1/2}$", L"$E_0$"))
axis((-3,2,0,4))
grid(true)
savefig("E_plots.pdf")

figure(2)
y = range(-1, 1, length=201)
x = ( 1 .+ y ) ./ ( 1 .- y )
plot(y, exp.(-x), y, erfcx.(x), y, 1 ./ ( 1 .+ x ))
legend((L"$g_1$", L"$g_{1/2}$", L"$g_0$"))
xlabel("y")
grid(true)
savefig("g_plots.pdf")
