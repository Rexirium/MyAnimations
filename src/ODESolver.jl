module ODESolver

export ODE_Problem, OnedimProblem, MultidimProblem, SecondOrderProblem
export ode_solve, quasi_newton

function quasi_newton(func::Function, xini::Tuple{T, T}; tol=1e-14, xtol=1e-12,  maxiter::Int=100) where T <: Real
    x0, x1 = xini
    i = 0
    f0, f1 = func(x0), func(x1)
    while func(x1) > tol || abs(x1 - x0) > xtol
        xnew = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, x1 = x1, xnew
        f0, f1 = f1, func(xnew)
        i += 1
        if i > maxiter
            error("Maximum iterations reached in quasi-Newton method.")
            break
        end
    end
    return x1
end

abstract type ODE_Problem end

struct OnedimProblem{T <: Number, R <: Real} <: ODE_Problem
    func::Function
    y0::T
    ts::LinRange{R}
    dt::R
end

struct MultidimProblem{T <: Number, R <: Real} <: ODE_Problem
    func::Function
    y0::Vector{T}
    ts::LinRange{R}
    dt::R
end

struct SecondOrderProblem{T <: Number, R <: Real} <: ODE_Problem
    func::Function
    y0::T
    dy0::T
    ts::LinRange{R}
    dt::R
end

function OnedimProblem(func::Function, y0::T, tspan::Tuple{R, R}, nsteps::Int=100) where {T <: Number, R <: Real}
    dt = (tspan[2] - tspan[1]) / (nsteps - 1)
    ts = LinRange{R}(tspan..., nsteps)
    OnedimProblem{T, R}(func, y0, ts, dt)
end
function OnedimProblem(func::Function, y0::T, tspan::Tuple{R, R}, dt::AbstractFloat) where {T <: Number, R <: Real}
    nsteps = floor(Int, (tspan[2] - tspan[1]) / dt) + 1
    ts = LinRange{R}(tspan..., nsteps)
    OnedimProblem{T, R}(func, y0, ts, dt)
end

function MultidimProblem(func::Function, y0::Vector{T}, tspan::Tuple{R, R}, nsteps::Int=100) where {T<:Number, R<:Real}
    dt = (tspan[2] - tspan[1]) / (nsteps - 1)
    ts = LinRange{R}(tspan..., nsteps)
    MultidimProblem{T, R}(func, y0, ts, dt)
end
function MultidimProblem(func::Function, y0::Vector{T}, tspan::Tuple{R, R}, dt::AbstractFloat) where {T<:Number, R<:Real}
    nsteps = floor(Int, (tspan[2] - tspan[1]) / dt) + 1
    ts = LinRange{R}(tspan..., nsteps)
    MultidimProblem(func, y0, ts, dt)
end

function SecondOrderProblem(func::Function, y0::T, dy0::T, tspan::Tuple{R,R}, nsteps::Int=100) where {T<:Number, R<:Real}
    dt = (tspan[2] - tspan[1]) / (nsteps - 1)
    ts = LinRange{R}(tspan..., nsteps)
    SecondOrderProblem{T, R}(func, y0, dy0, ts, dt)
end
function SecondOrderProblem(func::Function, y0::T, dy0::T, tspan::Tuple{R,R}, dt::AbstractFloat) where {T<:Number, R<:Real}
    nsteps = floor(Int, (tspan[2] - tspan[1]) / dt) + 1
    ts = LinRange{R}(tspan..., nsteps)
    SecondOrderProblem{T, R}(func, y0, dy0, ts, dt)
end


function ode_solve(prob::OnedimProblem{T, R}) where {T <: Number, R <: Real}
    h = prob.dt
    func = prob.func
    h_half = h / 2

    results = Vector{T}(undef, length(prob.ts))
    y = prob.y0

    ks = Vector{T}(undef, 4)
    for (i, t) in enumerate(prob.ts)
        results[i] = y

        ks[1] = func(t, y)
        ks[2] = func(t + h_half, y + h_half * ks[1])
        ks[3] = func(t + h_half, y + h_half * ks[2])
        ks[4] = func(t + h, y + h * ks[3])

        y += (h / 6) * (ks[1] + 2 * ks[2] + 2 * ks[3] + ks[4])
    end
    return results
end

function ode_solve(prob::MultidimProblem{T, R}) where {T <: Number, R <: Real}
    h = prob.dt
    funcs = prob.func
    h_half = h / 2

    results = Matrix{T}(undef, length(prob.y0), length(prob.ts))
    ys = copy(prob.y0)

    k1 = similar(ys)
    k2 = similar(ys)
    k3 = similar(ys)
    k4 = similar(ys)
    temp = similar(ys)

    for (i, t) in enumerate(prob.ts)
        results[:, i] .= ys

        k1 .= funcs(t, ys)
        temp = ys + h_half * k1
        k2 .= funcs(t + h_half, temp)
        temp = ys + h_half * k2
        k3 .= funcs(t + h_half, temp)
        temp = ys + h * k3
        k4 .= funcs(t + h, temp)

        ys += (h/6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end
    return results
end

function ode_solve(prob::SecondOrderProblem{T, R}) where {T <: Number, R <: Real}
    h = prob.dt
    func = prob.func
    h_half = h / 2

    nsteps = length(prob.ts)
    var_results = Vector{T}(undef, nsteps)
    deriv_results = Vector{T}(undef, nsteps)

    y, dy = prob.y0, prob.dy0

    ks = Vector{T}(undef, 4)
    kds = Vector{T}(undef, 4)
    for (i, t) in enumerate(prob.ts)
        var_results[i] = y
        deriv_results[i] = dy

        ks[1] = dy
        kds[1] = func(t, y, dy)
        ks[2] = dy + h_half * kds[1]
        kds[2] = func(t + h_half, y + h_half * ks[1], dy + h_half * kds[1])
        ks[3] = dy + h_half * kds[2]
        kds[3] = func(t + h_half, y + h_half * ks[2], dy + h_half * kds[2])
        ks[4] = dy + h * kds[3]
        kds[4] = func(t + h, y + h * ks[3], dy + h * kds[3])

        y += (h/6) * (ks[1] + 2*ks[2] + 2*ks[3] + ks[4])
        dy += (h/6) * (kds[1] + 2*kds[2] + 2*kds[3] + kds[4])
    end
    
    return var_results, deriv_results
end

end # module ODESolver