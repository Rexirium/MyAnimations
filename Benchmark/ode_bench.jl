using Revise
if !isdefined(Main, :ODESolver)
    include("../src/ODESolver.jl")
    using .ODESolver
end

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

    initial_state = [r0, θ0, dr0]
    tspan = (0.0, final_time)

    ts = range(tspan... , nsteps+1)
    prob = MultidimProblem(equations, initial_state, tspan, nsteps+1)
    @time results = ode_solve(prob)
    nothing
end