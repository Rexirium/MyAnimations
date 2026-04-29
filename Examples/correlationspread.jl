using CairoMakie

if !isdefined(Main, :Integral)
    include("../src/Integral.jl")
    using .Integral
end

let 
    J, Δ, μ = 1.0, 2.0, 1.0

    ens0(k) = -2J * cos(k) - μ
    ens(k) = sqrt(ens0(k)^2 + 4 * (Δ * sin(k))^2)
    ff(k::Real, x::Real) = 1/2 * cis(k*x) * (1.0 - ens0(k) / ens(k))
    f0(k::Real) = 1/2 * (1.0 - ens0(k) / ens(k))
    gg(k::Real, x::Real, t::Real) = -im * cis(- k*x - 2*ens0(k)*t) *(Δ * sin(k))/ens(k)

    nx = 500
    xs = range(-100, 100, nx)

    F2s = Vector{Float64}(undef, nx)
    for (i, x) in enumerate(xs)
        F2s[i] = abs2(integrate(ff, FilonCis(500, x), x; a=-π, b=π) / 2π)
    end
    F0 = abs2(integrate(f0, ChebyshevU(500); a=-π, b=π) / 2π)
    println("at x = 0 , |f(x)|² = $F0")

    fig = Figure()
    ax = Axis(fig[1,1], yscale=log10, 
        limits =(nothing, (1e-20, 1e0)))
    lines!(ax, xs, F2s)
    fig
end