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
