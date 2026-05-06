using CairoMakie

if !isdefined(Main, :Integral)
    include("../src/Integral.jl")
    using .Integral
end

let 
    J, Δ, μ = 1.0, 2.0, 1.0

    ens0(k) = -2J * cos(k) - μ
    ens(k) = sqrt(ens0(k)^2 + (2Δ * sin(k))^2)
    ff(k::Real, x::Real) = 1/2 * cis(k*x) * (1.0 - ens0(k) / ens(k))
    f0(k::Real) = 1/2 * (1.0 - ens0(k) / ens(k))
    gg(k::Real, x::Real, t::Real) = -im * cis(- k*x - 2*ens0(k)*t) *(Δ * sin(k))/ens(k)


    nx, nt = 500, 201
    
    xs = range(-100.0, 100.0, nx)
    ts = range(0.0, 20.0, nt)
    ts2 = range(0.1, 100, 1000)

    F2s = Vector{Float64}(undef, nx)
    G2s = Matrix{Float64}(undef, nx, nt)
    G2s2 = Vector{Float64}(undef, length(ts2))
    
    for (i, x) in enumerate(xs)
        F2s[i] = abs2(integrate(ff, FilonCis(1000, x), x; a=-π, b=π) / (2π))
    end
    #=
    for (j, t) in enumerate(ts)
        for (i, x) in enumerate(xs)
            G2s[i, j] = abs2(integrate(gg, FilonCis(500, x), -x, t; a=-π, b=π) / 2π)
        end
    end
    =#
    for (j, t) in enumerate(ts2)
        x = t
        G2s2[j] = abs2(integrate(gg, FilonCis(1000, x), -x, t; a=-π,b=π) / 2π)
    end
    
    set_theme!(Axis=(
        xtickalign=1,
        ytickalign=1, 
        xgridvisible=false,
        ygridvisible=false,
        xlabelsize=18,
        ylabelsize=18
    ), 
    Legend=(labelsize=18, framecolor=:grey, framewidth=1), 
    theme_latexfonts())
    #=
    fig = Figure()
    ax = Axis(fig[1,1], 
        xlabel=L"x", ylabel=L"t"
    )
    hm = heatmap!(ax, xs, ts, G2s, colorscale=sqrt)
    ablines!(ax, [0, 0], [1/4, -1/4], color=[:yellow, :yellow])

    Colorbar(fig[1,2], hm)
    fig
    save("Examples/figures/dynamics.png", fig)
    =#
    fig = Figure(size=(600, 700))
    ax1 = Axis(fig[1,1], yscale=log10,
        xlabel=L"x\ell" ,
        limits=(nothing, (1e-14, 1e0))   
    )
    lines!(ax1, xs, F2s, label=L"|f(\ell)|^2")
    axislegend(ax1, labelsize=18)

    ax2 = Axis(fig[2,1], yscale=log10,
        xlabel=L"t"    
    )
    lines!(ax2, ts2, G2s2, label=L"|g(\ell=vt, t)|^2")
    lines!(ax2, ts2, 1.2e-2 ./ ts2, label=L" \propto t^{-1}", linestyle=:dash)
    axislegend(ax2, labelsize=18)
    #save("Examples/figures/FandG.png", fig)
    fig
end