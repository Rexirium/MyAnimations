module ODESolver

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

function ode_solve(funcs::Function, y0::T, tspan::Tuple{R, R}; 
    nsteps::Int=100, h::Union{Real, Nothing}=nothing) where {T <: Real, R <: Real}
    
    if !isnothing(h)
        ti, tf, h = promote(tspan..., h)
        nsteps = ceil(Int, (tf - ti) / h)
    else
        ti, tf = tspan
        h = (tf - ti) / nsteps
    end
    h_half = h / 2

    results = Vector{T}(undef, nsteps + 1)
    t, y = ti, y0
    results[1] = y

    for i in 1:nsteps
        if t + h > tf
            h = tf - t
            h_half = h / 2
        end
        k1 = funcs(t, y)
        k2 = funcs(t + h_half, y + h_half * k1)
        k3 = funcs(t + h_half, y + h_half * k2)
        k4 = funcs(t + h, y + h * k3)

        y += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
        results[i + 1] = y
    end
    return results
end

function ode_solve(funcs::Function, y0::AbstractVector{T}, tspan::Tuple{R, R};
    nsteps::Int=100, h::Union{Real, Nothing}=nothing) where {T <: Number, R <: Real}
    if !isnothing(h)
        ti, tf, h = promote(tspan..., h)
        nsteps = ceil(Int, (tf - ti) / h)
    else
        ti, tf = tspan
        h = (tf - ti) / nsteps
    end
    h_half = h / 2

    results = Matrix{T}(undef, length(y0), nsteps + 1)
    ys = copy(y0)
    results[:, 1] .= ys
    t = ti

    k1 = similar(ys)
    k2 = similar(ys)
    k3 = similar(ys)
    k4 = similar(ys)
    temp = similar(ys)

    for i in 1:nsteps
        if t + h > tf
            h = tf - t
            h_half = h / 2
        end

        k1 .= funcs(t, ys)
        temp = ys + h_half * k1
        k2 .= funcs(t + h_half, temp)
        temp = ys + h_half * k2
        k3 .= funcs(t + h_half, temp)
        temp = ys + h * k3
        k4 .= funcs(t + h, temp)

        ys += (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        t += h
        results[:, i + 1] .= ys
    end
    return results
end

function ode_solve(func::Function, y0::T, dy0::T, tspan::Tuple{R, R}; 
    nsteps::Int=100, h::Union{Real, Nothing}=nothing) where {T <: Number, R <: Real}
    if !isnothing(h)
        ti, tf, h = promote(tspan..., h)
        nsteps = ceil(Int, (tf - ti) / h)
    else
        ti, tf = tspan
        h = (tf - ti) / nsteps
    end
    h_half = h / 2

    var_results = Vector{T}(undef, nsteps + 1)
    deriv_results = Vector{T}(undef, nsteps + 1)
    y, dy = y0, dy0
    var_results[1] = y
    deriv_results[1] = dy
    t = ti

    for i in 1:nsteps
        if t + h > tf
            h = tf - t
            h_half = h / 2
        end

        k1y = dy
        k1dy = func(t, y, dy)
        k2y = dy + h_half * k1dy
        k2dy = func(t + h_half, y + h_half * k1y, dy + h_half * k1dy)
        k3y = dy + h_half * k2dy
        k3dy = func(t + h_half, y + h_half * k2y, dy + h_half * k2dy)
        k4y = dy + h * k3dy
        k4dy = func(t + h, y + h * k3y, dy + h * k3dy)

        y += (h/6) * (k1y + 2*k2y + 2*k3y + k4y)
        dy += (h/6) * (k1dy + 2*k2dy + 2*k3dy + k4dy)
        t += h
        var_results[i + 1] = y
        deriv_results[i + 1] = dy
    end
    return var_results, deriv_results
end

end # module ODESolver