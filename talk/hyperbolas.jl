using PyPlot

r = 0.6
s = 0.7
ϕ = 0.8
μ = 1.0

@assert 0 < -s+ϕ < π/2
@assert 0 <  r+ϕ < π/2

u0 = range(-1.8, 1.8, length=50)
ur = range(-1.7, 1.7, length=50)
us = range(-1.5, 1.5, length=50)

figure(1)
x0 = μ * ( 1 .- cosh.(u0)*sin(ϕ) )
y0= μ * sinh.(u0) * cos(ϕ)
xr = μ * ( 1 .- cosh.(ur)*sin(r+ϕ) )
yr= μ * sinh.(ur) * cos(r+ϕ)
xs = μ * ( 1 .- cosh.(us)*sin(-s+ϕ) )
ys= μ * sinh.(us) * cos(-s+ϕ)
plot(xr, yr, x0, y0, xs, ys)
legend((L"$r$", L"$0$", L"$s$"))
grid(true)
axis("equal")
