function assess(α::T, β::T, x::T, tol::T) where T <: AbstractFloat
    Nopt = ceil(Integer, x^(1/α))
    term = NaN
    Nmin = -1
    for n = 1:Nopt
        term = exp( loggamma(n*α-β+1) - n*log(x) )
        if term < tol
            Nmin = n
            break
        end
    end
    return Nmin, term
end

