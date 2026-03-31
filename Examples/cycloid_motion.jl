using Revise
using GLMakie
if !isdefined(Main, :ODESolver)
    include("../src/ODE_Solver/ODESolver.jl")
    using .ODESolver
end

let 
    q, m = 1.0, 1.0
    E, B = 1.0, 1.0
    g = 0.02
    v0 = π

    ω = q * B / m
    γ = g / m
    ϵ = q * E / m

    T = 2π / ω
    R = (v0 + E / B) / ω

    function equation(t, u)
        du = similar(u)
        du[1] = u[3]
        du[2] = u[4]
        du[3] = ω * u[4] - γ * u[3]
        du[4] = -ω * u[3] - γ * u[4] + ϵ
        return du
    end

    final_time = 2.5T
    nsteps = 2000
    time_len = 10
    fps = 25

    initial = zeros(4)
    initial[3] = -v0
    tspan = (0.0, final_time)
    results = ode_solve(equation, initial, tspan; nsteps=nsteps)

    xs = results[1, :]
    ys = results[2, :]

    xmax = maximum(xs)

    xlo, xup = -R, xmax
    ylo, yup = -R, 3R

    nframes = fps * time_len
    step_per_frame = nsteps ÷ nframes

    set_theme!(theme_black())
    trajectory = Point2f[]
    current_loc = Point2f(initial[1:2])
    fig, ax, l = lines(trajectory, linewidth=2,
        axis = (; 
            aspect=DataAspect(),
            xlabel = L"x", xtickalign=1,
            ylabel = L"y", ytickalign=1,
            limits = (xlo, xup, ylo, yup))
    )
    s = scatter!(ax, [current_loc], color=:red, markersize=8)

    record(fig, "Examples/videos/cycloid_motion.mp4", 0:nframes; framerate=fps) do frame
        if frame > 0
            start = (frame - 1)*step_per_frame + 1
            stop = start + step_per_frame - 1

            new_points = Point2f.(xs[start:stop], ys[start:stop])
            append!(trajectory, new_points)
            Makie.update!(l, arg1 = trajectory) # Update the line with the new trajectory points
            Makie.update!(s, arg1 = [new_points[end]]) # Update the scatter point
            
        end
    end

end