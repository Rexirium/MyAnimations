using Revise
using GLMakie
include("../src/ODE_Solver/runge_kutta.jl")
include("attractors.jl")

let

    u0 = [1.0, 1.0, 1.0]

    final_time = 120.0
    h = 0.01
    time_len = 40
    fps = 50

    attractor1 = LorenzAttractor{Float64}(10.0, 28.0, 8/3)
    lorenz = get_ode_function(attractor1)

    tspan = (0.0, final_time)
    results = ode_solve(lorenz, u0, tspan; h=h)
    # ts = tspan[1]:h:tspan[2]

    nsteps = size(results, 2)
    nframes = time_len * fps
    step_per_frame = nsteps ÷ nframes

    points = Point3f[u0]
    current_point = Point3f(u0)
    colors = Int[0]

    set_theme!(theme_black())
    fig, ax, l = lines(points, color = colors, linewidth=1,
        colormap=:inferno, transparency=true,
        axis = (; type = Axis3, protrusions = (0, 0, 0, 0),
              viewmode = :fit, limits = (-30, 30, -30, 30, 0, 50))
    )
    s = scatter!(ax, [current_point], color=:white, markersize=8)

    record(fig, "Examples/videos/lorenz_attractor.mp4", 0:nframes; framerate=fps) do frame
        if frame > 0
            start = (frame - 1)*step_per_frame + 1
            stop = start + step_per_frame - 1

            new_points = Point3f.(eachcol(results[:, start : stop]))
            append!(points, new_points)
            append!(colors, fill(frame, step_per_frame))

            ax.azimuth[] = 1.7π + 0.5 * sin(2π * frame / nframes)
            Makie.update!(l, arg1 = points, color = colors) # Makie 0.24+
            Makie.update!(s, arg1 = [new_points[end]]) # Update the scatter point to the current position
            l.colorrange = (0, frame)
        end
    end
end