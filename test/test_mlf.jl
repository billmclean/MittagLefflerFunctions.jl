
include("mlf.jl")

args = [ (1.0, 1.0, 1.0),
         (0.5, 2.0, -0.2),
         (1.5, 1.5, 2.0),
         (2.5, 1.6, 3.0),
         (0.7, 1.1, 2.4) ]

@printf("%8s  %8s  %8s  %12s\n\n", "α", "β", "x", "E_αβ(x)")
for arg_list in args
    α, β, x = arg_list
    y = mlf(α, β, x, 10)
    @printf("%8.4f  %8.4f  %8.4f  %12.8f\n", α, β, x, y)
end
