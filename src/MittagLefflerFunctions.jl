module MittagLefflerFunctions

import SpecialFunctions

Γ = SpecialFunctions.gamma

"""
    rΓ(a)

Returns the reciprocal Gamma function `1/Γ(a)` for any real `a`.
"""
function rΓ(a::T) where T <: AbstractFloat
    if a > 0
        return one(T) / Γ(a)
    else
        return sinpi(a) * Γ(1-a) / π
    end
end

include("submodules/Pade.jl")
include("submodules/Quadrature.jl")
include("submodules/MLF.jl")

using .Pade
using .Quadrature
using .MLF

export MLPade
export MLQuad
export mlf

end # module
