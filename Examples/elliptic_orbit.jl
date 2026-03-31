using Revise
using GLMakie
include("../src/ODE_Solver/runge_kutta.jl")

let 
    m, M = 1.0, 10.0
    Ω = 4.0
    λ = 4.0
    r0, v0 = 1.0, 2.0
    θ0 = 0.0
    dr0 = 0.0
    E = 0.5* m * v0*v0 - M * m/ r0
    L = m * r0 * v0
    a = - M * m / (2*E)
    T = 2π / Ω
    x0, y0 = r0 * cos(θ0), r0 * sin(θ0)

    Vp(r) =  m * (Ω * r)^2 / 2 * exp( - (r / λ)^2)
    Fp(r) = - m * Ω * Ω * r * (1.0 - (r/λ)^2) * exp( - (r/λ)^2)
    function equations(t, state)
        r = state[1]
        dr = state[3]
        d2r = L*L/(m*m*r^3) + Fp(r)/m
        dθ = L/(m*r*r)
        return [dr, dθ, d2r]
    end

    final_time = 5T
    nsteps = 5000
    time_len = 20
    fps = 25

    initial_state = [r0, θ0, dr0]
    tspan = (0.0, final_time)

    ts = range(tspan... , nsteps+1)

    results = ode_solve(equations, initial_state, tspan; nsteps=nsteps)

    xs = results[1, :] .* cos.(results[2, :])
    ys = results[1, :] .* sin.(results[2, :])

    xmin, xmax = extrema(xs)
    ymin, ymax = extrema(ys)

    xcenter = (xmin + xmax) / 2
    ycenter = (ymin + ymax) / 2
    deltax, deltay = xmax - xmin, ymax - ymin
    margin = 1.5
    xlo, xup = xcenter - margin * deltax / 2, xcenter + margin * deltax / 2
    ylo, yup = ycenter - margin * deltay / 2, ycenter + margin * deltay / 2 

    nframes = fps * time_len
    step_per_frame = nsteps ÷ nframes

    set_theme!()
    trajectory = Observable(Point2f[])
    solar_loc = Point2f(0.0, 0.0)
    planet_loc = Observable(Point2f(x0, y0))
    radius = Observable(Point2f[(0, 0), (x0, y0)])
    time_text = Observable("time = 0.00 T")

    fig = Figure()
    ax = Axis(fig[1, 1], title=time_text, 
        titlesize=18,
        aspect=DataAspect(),
        xlabel=L"x",
        ylabel=L"y",
        limits=(xlo, xup, ylo, yup),
        xtickalign=1,
        ytickalign=1)

    lines!(ax, trajectory, color=:purple, linewidth=1, label="orbit")
    lines!(ax, radius, color=:gray, linewidth=1, linestyle=:dash)
    scatter!(ax, solar_loc, markersize=25, color=:orange, label="Sun")
    scatter!(ax, planet_loc, markersize=10, color=:blue, label="Planet")
    
    axislegend()

    new_traj = Vector{Point2f}(undef, step_per_frame)
    record(fig, "Examples/videos/elliptic_orbit.mp4", 0:nframes; framerate=fps) do frame
        if frame > 0
            start = (frame - 1)*step_per_frame + 1
            stop = start + step_per_frame - 1
            idx = stop + 1

            @. new_traj = Point2f(xs[start:stop], ys[start:stop])

            planet_loc[] = Point2f(xs[idx], ys[idx])
            radius[] = Point2f[(0.0, 0.0), (xs[idx], ys[idx])]
            trajectory[] = append!(trajectory[], new_traj)
            time_text[] = "time = $(rpad(round(ts[idx]/T; digits = 2), 4, "0")) T"
        end
    end

end