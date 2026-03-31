abstract type AbstractAttractor end

mutable struct LorenzAttractor{T <: Real} <: AbstractAttractor
    σ::T
    ρ::T
    β::T
end

mutable struct RosslerAttractor{T <: Real} <: AbstractAttractor
    a::T
    b::T
    c::T
end

mutable struct ChenAttractor{T <: Real} <: AbstractAttractor
    a::T
    b::T
    c::T
end

mutable struct AizawaAttractor{T <: Real} <: AbstractAttractor
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
end

function get_ode_function(attr::LorenzAttractor)
    function lorenz(t, u)
        du = similar(u)
        du[1] = attr.σ * (u[2] - u[1])
        du[2] = u[1] * (attr.ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - attr.β * u[3]
        return du
    end
    return lorenz
end

function get_ode_function(attr::RosslerAttractor)
    function rossler(t, u)
        du = similar(u)
        du[1] = -u[2] - u[3]
        du[2] = u[1] + attr.a * u[2]
        du[3] = attr.b + u[3] * (u[1] - attr.c)
        return du
    end
    return rossler
end

function get_ode_function(attr::ChenAttractor)
    function chen(t, u)
        du = similar(u)
        du[1] = attr.a * (u[2] - u[1])
        du[2] = (attr.c - attr.a) * u[1] - u[1] * u[3] + attr.c * u[2]
        du[3] = u[1] * u[2] - attr.b * u[3]
        return du
    end
    return chen
end

function get_ode_function(attr::AizawaAttractor)
    function aizawa(t, u)
        du = similar(u)
        du[1] = (u[3] - attr.b) * u[1] - attr.d * u[2]
        du[2] = attr.d * u[1] + (u[3] - attr.b) * u[2]
        du[3] = attr.c + attr.a * u[3] - (u[3]^3) / 3 
            - (u[1] * u[1] + u[2] * u[2]) * (1 + attr.e * u[3]) + attr.f * u[3] * (u[1]^3)
        return du
    end
    return aizawa
end