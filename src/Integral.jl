module Integral

using FastGaussQuadrature: gausslegendre, gausshermite, gausslaguerre, gaussjacobi

export Chebyshev, ChebyshevT, ChebyshevU, ChebyshevV, ChebyshevW
export LegendreP, HermiteH, LaguerreL, JacobiP
export Filon, FilonCis, FilonCos, FilonSin
export Trapezoid, SimpsonRule, integrate

abstract type AbstractQuadrature end

abstract type Chebyshev <: AbstractQuadrature end
abstract type Filon <: AbstractQuadrature end

struct Trapezoid <: AbstractQuadrature num::Int end
struct SimpsonRule <: AbstractQuadrature num::Int end

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

function integrate(integrand::Function, trap::Trapezoid, args...; a::Real=-1.0, b::Real=1.0)
    h = (b - a) / trap.num
    total = 0.5 * (integrand(a, args...) + integrand(b, args...))
    xi = a + h
    for _ in 1 : trap.num - 1
        total += integrand(xi, args...)
        xi += h
    end
    return h * total
end

function integrate(integrand::Function, simp::SimpsonRule, args...; a::Real=-1.0, b::Real=1.0)
    n = 2simp.num
    h = (b - a) / n
    total = integrand(a, args...) + integrand(b, args...)
    xi = a + h
    for i in 1 : n - 1
        wt = isodd(i) ? 4 : 2
        total += wt * integrand(xi, args...)
        xi += h
    end
    return h/3 * total
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
    θ3 = θ^3
    θ2 = θ^2
    if abs(θ) > 1/6
        alpha = (θ*θ + θ/2 * sin(2θ) + cos(2θ) - 1) / θ3
        beta = (θ*(3 + cos(2θ)) - 2 * sin(2θ)) / θ3
        gamma = 4 * (sin(θ) - θ*cos(θ)) / θ3
    else
        alpha = θ3 * evalpoly(θ2, (2/45, - 2/315, 2/4725))
        beta = evalpoly(θ2, (2/3, 2/15, -4/105, 2/567))
        gamma = evalpoly(θ2, (4/3, -2/15, 1/210, -1/11340))
    end
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


