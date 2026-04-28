abstract type IntegralAlgorithm end

struct ChebyshevT <: IntegralAlgorithm
    order::Int
end


function integrate(integrand::Function, alg::ChebyshevT, args...; xmin::Real=-1.0, xmax::Real=1.0)
    mid = (xmin + xmax) / 2
    half = (xmax - xmin) / 2
    order = alg.order
    weight = π / order

    total = zero(weight)
    for i in 1:order
        root = cospi((2i - 1) / (2order))
        wfun = sinpi((2i - 1) / (2order))
        xi = mid + half * root

        total += weight * integrand(xi, args...) * wfun
    end
    return half * total
end

let 
    ff(x, a, b) = exp(-a*x) * sin(b*x)
    @time integrate(ff, ChebyshevT(10), 1.0, 2.0; xmin=0.0, xmax=2π)
end