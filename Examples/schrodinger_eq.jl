using Revise
using GLMakie
if !isdefined(Main, :PDESolver)
    include("../src/PDE_Solver/PDESolver.jl")
    using .PDESolver
end

let 
    m = 1.0
    xa, xb = -1.0, 1.0
    tf = 0.05
    nx, nt = 1000, 48000

    sigma = 0.02
    x0, p0 = -0.5, 200.0

    function initial(x)
        normfactor = 1/(2π)^(1/4) / sqrt(sigma)
        phase = complex(cos(p0 * x), sin(p0 * x))
        return normfactor * exp( - 1/4 * (x - x0)*(x - x0)/sigma/sigma) * phase
    end
    potential(x) = 2*p0*p0/m * x*x

    prob = StaticSchrodingerProblem{Float64}(m, potential, initial, (xa, xb), tf)
    bc = DirichletCondition(0.0, 0.0)
    sol = pde_solve(prob, bc, nx, nt)

    # Visualization
    fps = 24
    nframes = 20 * fps
    step_per_frame = nt ÷ nframes

    step = Observable(1)
    xs = range(xa, xb, nx+1)
    ts = range(0.0, 1e3 * tf, nt+1)
    rho = @lift begin
       psi = sol[:, $step]
       abs2.(psi) 
    end

    fig = lines(xs, rho, linewidth=2,
        axis = (title=@lift("t = $(rpad(round(ts[$step], digits=2), 4, "0")) ms"),
                limits=(nothing, (-2, 25))))

    record(fig, "Examples/videos/wavepack.mp4", 0:nframes) do frame
        step[] = step_per_frame * frame + 1
    end    
    
end