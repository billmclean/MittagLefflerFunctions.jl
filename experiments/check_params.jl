
function error_terms(N)
    ϕ = 1.172104
    a = 1.081792
    b = 4.492075
    h = a / N
    μ = b * N
    E1 = exp(-π*(π-2ϕ)/h)
    E2 = exp(μ-2π*ϕ/h)
    E3 = exp(μ*(1-cosh(a)*sin(ϕ)))
    return E1, E2, E3
end


