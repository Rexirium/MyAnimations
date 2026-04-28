using CairoMakie

if !isdefined(Main, :Integral)
    include("../src/Integral.jl")
    using .Integral
end

let 
    J, Δ, μ = 1.0, 2.0, 1.0

    eps0(k) = -2J * cos(k) - μ
    eps(k) = sqrt(eps0(k)^2 + 4 * (Δ * sin(k))^2)
    ff(k::Real, x::Real) = 1/2 * cis(k*x) * (1.0 - eps0(k) / eps(k))
    gg(k::Real, x::Real, t::Real) = -im * cis(- k*x - 2*eps0(k)*t) *(Δ * sin(k))/eps(k)

    nx = 1001
    xs = range(-100, 100, nx)

    F2s = Vector{Float64}(undef, nx)
    for (i, x) in enumerate(xs)
        F2s[i] = abs2(integrate(gg, FilonCis(500, -x), -x, 20.0; a=-π, b=π) / 2π)
    end
    lines(xs, F2s)
end