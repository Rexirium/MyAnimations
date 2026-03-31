using GLMakie

include("../src/PDE_Solver/convection_diffusion.jl")

let 
    xa, xb = -1.0, 1.0
    tf = 5e-2
    nx, nt = 500, 24000
    
    diff(x) = one(x)
    x0, σ = -0.5, 0.02
    h = 1 / sqrt(2π) / σ 
    initial(x) = h * exp( - (x-x0)*(x-x0) / (2 * σ * σ))

    prob = ConstConvDiffProblem(50.0, 1.0, 100.0, initial, (xa, xb), tf)
    bc = DirichletCondition(0.0, 0.0)

    sol = pde_solve(prob, bc, nx, nt)

    # Visualization
    fps = 24
    nframes = 20 * fps
    step_per_frame = nt ÷ nframes

    step = Observable(1)
    xs = range(xa, xb, nx+1)
    ts = range(0.0, 1e3 * tf, nt+1)
    us = @lift(sol[:, $step])

    fig = lines(xs, us, linewidth=2,
        axis = (limits=(nothing, (-0.2*h, 1.2*h)), 
            title = @lift("t = $(rpad(round(ts[ $step ], digits=2), 4, "0")) ms"),))

    record(fig, "Examples/videos/gaussian_diffusion.mp4", 0:nframes) do frame
        step[] = step_per_frame * frame + 1
    end    
end