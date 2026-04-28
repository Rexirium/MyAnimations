module Integral

using FastGaussQuadrature: gausslegendre, gausshermite, gausslaguerre, gaussjacobi

export Chebyshev, ChebyshevT, ChebyshevU, ChebyshevV, ChebyshevW
export LegendreP, HermiteH, LaguerreL, JacobiP
export Filon, FilonCis, FilonCos, FilonSin
export integrate

abstract type AbstractQuadrature end
abstract type Chebyshev <: AbstractQuadrature end
abstract type Filon <: AbstractQuadrature end

struct ChebyshevT <: Chebyshev order::Int end
struct ChebyshevU <: Chebyshev order::Int end
struct ChebyshevV <: Chebyshev order::Int end
struct ChebyshevW <: Chebyshev order::Int end

struct LegendreP <: AbstractQuadrature order::Int end
struct HermiteH <: AbstractQuadrature order::Int end
struct LaguerreL <: AbstractQuadrature order::Int end

struct JacobiP <: AbstractQuadrature
    order::Int
    alpha::Float64
    beta::Float64
end

struct FilonCis <: Filon 
    num::Int
    freq::Float64
end
struct FilonCos <: Filon
    num::Int
    freq::Float64
end
struct FilonSin <: Filon
    num::Int
    freq::Float64
end

function chebyshev(cheby::ChebyshevT, i::Int)
    n = cheby.order
    return sincospi((2i - 1)/(2n))..., π / n
end

function chebyshev(cheby::ChebyshevU, i::Int)
    n1 = cheby.order + 1
    wfi, xi = sincospi(i / n1)
    return 1/wfi, xi, π * (1 - xi*xi) / n1
end

function chebyshev(cheby::ChebyshevV, i::Int)
    nodd = 2cheby.order + 1
    xi = cospi((2i - 1) / nodd)
    wfi = tanpi((i - 0.5) / nodd)
    return wfi, xi, 2π * (1 + xi) / nodd
end

function chebyshev(cheby::ChebyshevW, i::Int)
    nodd = 2cheby.order + 1
    xi = cospi((2i) / nodd)
    wfi = cot(π * i / nodd)
    return wfi, xi, 2π * (1 - xi) / nodd
end

function integrate(integrand::Function, cheby::Chebyshev, args...; a::Real=-1.0, b::Real=1.0)
    mid = (a + b) / 2
    half = (b - a) / 2

    total = zero(integrand(mid, args...))
    for i in 1 : cheby.order
        wfun, root, weight = chebyshev(cheby, i)
        xi = mid + half * root

        total += weight * integrand(xi, args...) * wfun
    end
    return half * total
end

function integrate(integrand::Function, chebyt::ChebyshevT, args...; a::Real=-1.0, b::Real=1.0)
    mid = (a + b) / 2
    half = (b - a) / 2
    n = chebyt.order
    weight = π / n

    total = zero(integrand(mid, args...))
    for i in 1 : n
        wfun, root = sincospi((i - 0.5) / n)
        xi = mid + half * root

        total += weight * integrand(xi, args...) * wfun
    end
    return half * total
end

function integrate(integrand::Function, leg::LegendreP, args...; a::Real=-1.0, b::Real=1.0)
    mid = (a + b) / 2
    half = (b - a) / 2

    n = leg.order
    roots, weights = gausslegendre(n)

    total = zero(integrand(mid, args...))
    for i in 1 : n
        xi = mid + half * roots[i]
        total += weights[i] * integrand(xi, args...)
    end
    return half * total
end

function integrate(integrand::Function, herm::HermiteH, args...)
    n = herm.order
    roots, weights = gausshermite(n)

    total = zero(integrand(0.0, args...))
    for i in 1 : n
        xi = roots[i]
        total += weights[i] * integrand(xi, args...) * exp(xi * xi)
    end
    return total
end

function integrate(integrand::Function, lag::LaguerreL, args...; a::Real=0.0)
    n = lag.order
    roots, weights = gausslaguerre(n)

    total = zero(integrand(1.0, args...))
    for i in 1 : n
        xi = roots[i]
        total += weights[i] * integrand(xi + a, args...) * exp(xi) 
    end
    return total
end

function integrate(integrand::Function, jcb::JacobiP, args...; a::Real=-1.0, b::Real=1.0)
    mid = (a + b) / 2
    half = (b - a) / 2

    n = jcb.order
    a, b = jcb.alpha, jcb.beta
    roots, weights = gaussjacobi(n, a, b)

    total = zero(integrand(mid, args...))
    for i in 1 : n
        ri = roots[i]
        xi = mid + half * ri

        wfi = (1.0 - ri)^(-a) * (1.0 + ri)^(-b)
        total += weights[i] * integrand(xi, args...) * wfi
    end
    return half * total
end

function make_filon_coeff(θ::Float64)
    alpha = (θ*θ + θ/2 * sin(2θ) + cos(2θ) - 1) / θ^3
    beta = (θ*(3 + cos(2θ)) - 2 * sin(2θ)) / θ^3
    gamma = 4 * (sin(θ) - θ*cos(θ)) / θ^3
    return alpha, beta, gamma
end

endpoint(filon::FilonCis, fa, fb, a::Real, b::Real) = im * (fa - fb)
endpoint(filon::FilonCos, fa, fb, a::Real, b::Real) = begin
    k = filon.freq
    fb * tan(k * b) - fa * tan(k * a)
end
endpoint(filon::FilonSin, fa, fb, a::Real, b::Real) = begin
    k = filon.freq
    fa * cot(k * a) - fb * cot(k * b)
end

function integrate(integrand::Function, filon::Filon, args...; a::Real=-1.0, b::Real=1.0)
    n, k = filon.num, filon.freq
    h = (b - a) / (2n)
    θ = k * h
    α, β, γ = make_filon_coeff(θ)

    fa, fb = integrand(a, args...), integrand(b, args...)
    S_even = 0.5 * (fa + fb)

    S_odd = zero(S_even)
    xi = a + h
    for i in 1 : 2n - 1
        if isodd(i)
            S_odd += integrand(xi, args...)
        else
            S_even += integrand(xi, args...)
        end
        xi += h
    end
    val = endpoint(filon, fa, fb, a, b)
    
    return h * (α * val + β * S_even + γ * S_odd)
end

end # Integral


