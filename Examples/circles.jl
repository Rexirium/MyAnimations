using LinearAlgebra
using GLMakie

let 
    R, r = 2, 1
    Ω = 1.0
    ω = (R - r) / r * Ω
    k = gcd(R - r, r)
    T = 2π / Ω
    tf = (R - r) / k * T

    nsteps = 1000
    time_len = 10
    fps = 25

    ts = range(0.0, tf, nsteps + 1)

    xc = (R - r) * cos.(Ω * ts)
    yc = (R - r) * sin.(Ω * ts)
    
    xs = xc + r * cos.(ω * ts)
    ys = yc - r * sin.(ω * ts)

    vx = - (R - r) * Ω * sin.(Ω * ts) - r * ω * sin.(ω * ts)
    vy = (R - r) * Ω * cos.(Ω * ts) - r * ω * cos.(ω * ts)

    nframes = fps * time_len
    spf = nsteps ÷ nframes

    set_theme!(Axis=(
        titlesize=20, 
        xtickalign=1, ytickalign=1,
    ))
    idx = Observable(1)

    θ0 = π/6 : π/6 : 2π
    
    trajectory = @lift(Point2f.(xs[1 : $idx], ys[1 : $idx]))
    center_loc = @lift(Point2f(xc[$idx], yc[$idx]))
    radius = @lift(Point2f[(xc[$idx], yc[$idx]), (xs[$idx], ys[$idx])])
    point_loc = @lift(Point2f(xs[$idx], ys[$idx]))
    velocity_vec = @lift(Point2f(0.5 * vx[$idx], 0.5 * vy[$idx]))
    dotlist = @lift begin
        list = Point2f[]
        for θ in θ0
            θs = ω * ts[$idx] + θ
            push!(list, Point2f(xc[$idx] + r * cos(θs), yc[$idx] - r * sin(θs)))
        end
        list
    end

    fig = Figure()
    ax = Axis(fig[1, 1], 
        title="Hypotrochoid", 
        xlabel="x", ylabel="y",
        aspect=DataAspect()
    )
    arc!(ax, Point2f(0, 0), R, 0, 2π, linewidth=3, color=:black)
    arc!(ax, Point2f(0, 0), R - r, 0, 2π, linewidth=1, color=:gray, linestyle=:dash)

    arc!(ax, center_loc, r, 0, 2π, linewidth=2, color=:blue)
    lines!(ax, trajectory, color=:purple, linewidth=1.5)
    lines!(ax, radius, color=:black, linewidth=1)
    scatter!(ax, center_loc, color=:blue, markersize=12)
    scatter!(ax, dotlist, color=:blue, markersize=12)
    scatter!(ax, point_loc, color=:red, markersize=16)
    arrows2d!(ax, point_loc, velocity_vec, color=:green)

    new_traj = Vector{Point2f}(undef, spf)
    record(fig, "Examples/videos/hypotrochoid.gif", 0:nframes; framerate=fps) do i
        idx[] = min(i * spf + 1, nsteps + 1)
    end

end