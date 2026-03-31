using Revise
using GLMakie
include("../src/ODE_Solver/runge_kutta.jl")

let
    g = 9.8
    l1, l2 = 1.0, 1.0
    m1, m2 = 1.0, 1.0

    θ1_0, ω1_0 = 0.0, 0.0 
    θ2_0, ω2_0 = 0.0, 8.0

    final_time = 10.0
    h = 0.01
    time_len = 10
    fps = 25

    tspan = (0.0, final_time)
    initial_state = [θ1_0, ω1_0, θ2_0, ω2_0]

    function equations(t, state)
        θ1, ω1, θ2, ω2 = state
        dθ1, dθ2 = ω1, ω2

        delta = θ1 - θ2
        denom = 2m1 + m2 - m2*cos(2*delta)
        dω1 = (- (2m1 + m2)*g*sin(θ1) - m2*g*sin(θ1-2θ2)
            - 2*m2*sin(delta)*(ω2^2*l2 + ω1^2*l1*cos(delta))) / (l1 * denom)
        dω2 = (2*sin(delta)*(ω1^2*l1*(m1 + m2) + g*(m1 + m2)*cos(θ1)
            + ω2^2*l2*m2*cos(delta))) / (l2 * denom)

        return [dθ1, dω1, dθ2, dω2]
    end

    results = ode_solve(equations, initial_state, tspan; h=h)
    ts = tspan[1] : h : tspan[2]
    if ts[end] < tspan[2]
        push!(ts, tspan[2])
    end

    x1 = l1 * sin.(results[1, :])
    y1 = -l1 * cos.(results[1, :])
    x2 = x1 .+ l2 * sin.(results[3, :])
    y2 = y1 .- l2 * cos.(results[3, :])

    nsteps = size(results, 2)
    nframes = time_len * fps
    step_per_frame = nsteps ÷ nframes

    trajectory = Observable(Point2f[])
    pend_lines = Observable(Point2f[(0, 0), (0, -1), (0, -2)])
    pend_points = Observable(Point2f[(0, 0), (0, -1), (0, -2)])
    time_texts = Observable("Double pendulum at time = 0.00 s")

    fig = Figure()
    ax = Axis(fig[1, 1], title=time_texts, 
        aspect=DataAspect(),
        xlabel="x", 
        ylabel="y", 
        limits=(-2.5, 2.5, -2.5, 1.5), 
        xticks=(-2:2:2),
        yticks=(-2:2:2),
        xtickalign=1, 
        ytickalign=1 
        )
    

    lines!(ax, trajectory, color=:orange, linewidth=1)
    lines!(ax, pend_lines, linewidth=3)
    scatter!(ax, pend_points, color=:red, markersize=15)

    new_traj = Vector{Point2f}(undef, step_per_frame)
    record(fig, "Examples/videos/double_pendulum.mp4", 0:nframes; framerate=fps) do frame

        if frame > 0 
            start = (frame - 1)*step_per_frame + 1
            stop = start + step_per_frame - 1
            idx = stop + 1

            current_points = Point2f.([0.0, x1[idx], x2[idx]], 
                [0.0, y1[idx], y2[idx]])
            
            pend_lines[] = current_points
            pend_points[] = current_points

            @. new_traj = Point2f(x2[start:stop], y2[start:stop])
            trajectory[] = append!(trajectory[], new_traj)

            time_texts[] = "Double pendulum at time = $(rpad(round(ts[idx]; digits=2), 4, "0")) s"
        end
    end

end