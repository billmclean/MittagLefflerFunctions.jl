using PyPlot
using SpecialFunctions

function g(y)
    x = (1+y) / (1-y)
    return ( erfcx(x) - one(y) ) / x
end

figure(1)
x = range(0, 5, length=201)
plot(x, erfcx.(x))

figure(2)
y = range(-1, 1, length=201)
plot(y, g.(y))
 
