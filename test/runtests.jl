using MittagLefflerFunctions
using Test
using SpecialFunctions: erfcx

@testset "mlf" begin
    for k = 1:3
        include("mlf/test$k.jl")
    end
end
