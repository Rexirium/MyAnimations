using Revise
if !isdefined(Main, :ODESolver)
    include("../src/ODE_Solver/ODESolver.jl")
    using .ODESolver
end
if !isdefined(Main, :Integral)
    include("../src/Integral.jl")
    using .Integral
end
using CairoMakie

elliptic(x::Real, k) = sqrt(1.0 - k^2 * sin(x)^2)

let 
    v0, R = 1.0, 1.0
    ff(t::Real, x, dx) = - 2 * (R * v0)^2 * x/(R*R + 4*x*x)^2
    tf = 20.0
    num = 500

    prob = SecondOrderProblem(ff, 0.0, -v0/2, (0.0, tf), num + 1)
    @time res = ode_solve(prob)

    ts = LinRange(0, tf, num + 1)
    θs = LinRange(0, 5π, num + 1)
    xs_anal = R/2 * sin.(θs)
    
    E_sqrt2 = integrate(elliptic, ChebyshevU(200), 1/√2 ; a=0, b=π/2)
    ts_anal = similar(ts)
    for (i, θ) in enumerate(θs)
        ts_anal[i] = √2 * R/v0 * (E_sqrt2 - integrate(elliptic, ChebyshevU(200), 1/√2 ; a=0, b=π/2 - θ))
    end
    
    
    set_theme!(Axis=(
        xtickalign=1, ytickalign=1,
        xlabelsize=18, ylabelsize=18, 
        xticklabelsize=16, yticklabelsize=16,
    ), theme_latexfonts())
    
    fig = Figure()
    ax = Axis(fig[1,1],
        xlabel=L"t", ylabel=L"x"
    )
    #lines!(ax, ts, res[1], label="tube oscillate")
    lines!(ax, ts,  res[1] .+ v0/2 * ts, label="ODE result")
    lines!(ax, ts_anal, - xs_anal .+ v0/2 * ts_anal, label="analytical result", linestyle=:dash)
    axislegend(ax)
    fig
    
end