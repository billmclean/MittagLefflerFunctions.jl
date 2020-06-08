module MittagLefflerFunctions

import SpecialFunctions

Γ = SpecialFunctions.gamma

include("submodules/Pade.jl")
include("submodules/Quadrature.jl")

using .Pade
using .Quadrature

export MLPade
export MLQuad1, MLQuad2, MLQuad3

end # module
