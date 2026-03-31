mutable struct ConvectionProblem{C<:Union{Fct, Real}, T<:Real} <: AbstractProblem
    velocity::C
    initial::Fct
    xlim::Tuple{T, T}
    tfin::T

    ConvectionProblem(velocity::Fct, initial::Fct, xlim::Tuple{T, T}, 
        tfin::T) where T<:Real = new{Fct, T}(velocity, initial, xlim, tfin)
    ConvectionProblem(velocity::R, initial::Fct, xlim::Tuple{T, T}, 
        tfin::T) where {R<:Real, T<:Real} = new{R, T}(velocity, initial, xlim, tfin)
end

mutable struct DiffusionProblem{D<:Union{Fct, Real}, T<:Real} <: AbstractProblem
    diffusivity::D
    initial::Fct
    xlim::Tuple{T, T}
    tfin::T

    DiffusionProblem(diffusivity::Fct, initial::Fct, xlim::Tuple{T, T}, 
        tfin::T) where T<:Real = new{Fct, T}(diffusivity, initial, xlim, tfin)
    DiffusionProblem(diffusivity::R, initial::Fct, xlim::Tuple{T, T}, 
        tfin::T) where {R<:Real, T<:Real} = new{R, T}(diffusivity, initial, xlim, tfin)
end

mutable struct ConstConvDiffProblem{C<:Real, T<:Real} <: AbstractProblem
    velocity::C
    diffusivity::C
    sourcerate::C
    initial::Fct
    xlim::Tuple{T, T}
    tfin::T

    ConstConvDiffProblem(velocity::C, diffusivity::C, sourcerate::C, initial::Fct, 
        xlim::Tuple{T, T}, tfin::T) where {C<:Real, T<:Real} = new{C, T}(
            velocity, diffusivity, sourcerate, initial, xlim, tfin
        )
end

mutable struct VariaConvDiffProblem{T<:Real} <: AbstractProblem
    velocity::Fct
    diffusivity::Fct
    sourcerate::Fct
    initial::Fct
    xlim::Tuple{T, T}
    tfin::T

    VariaConvDiffProblem(velocity::Fct, diffusivity::Fct, sourcerate::Fct, initial::Fct, 
        xlim::Tuple{T, T}, tfin::T) where T<:Real = new{T}(
            velocity, diffusivity, sourcerate, initial, xlim, tfin
        )
end

mutable struct TimedDiffusionProblem{T<:Real} <: AbstractProblem
    diffusivity::Fct
    initial::Fct
    xlim::Tuple{T, T}
    tfin::T

    TimedDiffusionProblem(diffusivity::Fct, initial::Fct, xlim::Tuple{T, T}, 
        tfin::T) where T<:Real = new{T}(diffusivity, initial, xlim, tfin)
end

function pde_solve(prob::ConvectionProblem{R, T}, bc::ConstantCondition, 
        nx::Int, nt::Int) where {R<:Real, T<:Real}
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cconv = dt / (4 * dx) * prob.velocity
    println("Courant number is: $(abs(4 * cconv))")
    if abs(4 * cconv) > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    us = prob.initial.(xs)

    dl = fill(cconv, nx-2)
    matrixret = Tridiagonal(-dl, ones(typeof(cconv), nx-1), dl)
    bs = Vector{eltype(us)}(undef, nx-1)

    solution = Matrix{eltype(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us

    for i in 2:nt+1
        bs .= us[2:nx] + cconv * (us[1:nx-1] - us[3:end])
        bs[1] += cconv * us[1]
        bs[end] -= cconv * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution
end

function convectionmatrix(c::Vector{T}) where T<:Real
    n = length(c)
    dd = -(c[3:n] - c[1:n-2])
    dl = c[3:n-1]
    du = - c[2:n-2]
    return Tridiagonal(dl, dd, du), c[2], - c[n-1]
end

function pde_solve(prob::ConvectionProblem{F, T}, bc::ConstantCondition, 
        nx::Int, nt::Int) where {F<:Fct, T<:Real}
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cc = dt / (4 * dx)
    println("Courant number is: $(4cc)")
    if 4cc > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    us = prob.initial.(xs)
    cconv = cc * prob.velocity.(xs)

    matrix, edgel, edgeu = convectionmatrix(cconv)
    matrixret = I(nx-1) .- matrix
    matrixadv = I(nx-1) .+ matrix
    bs = Vector{eltype(us)}(undef, nx-1)

    solution = Matrix{eltype(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us
    for i in 2:nt+1
        bs .= matrixadv * us[2:nx]
        bs[1] += 2 * edgel * us[1]
        bs[end] += 2 * edgeu * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution

end

function diffusionmatrix(c::Vector{T}) where T <: Real
    n = length(c) - 1
    dd = Vector{T}(undef, n-1)
    dl = Vector{T}(undef, n-1)
    du = Vector{T}(undef, n-1)
    for j in 1:n-1
        dd[j] = -2 * c[j+1]
        dl[j] = c[j+1] - (c[j+2] - c[j]) / 4
        du[j] = c[j+1] + (c[j+2] - c[j]) / 4
    end
    return Tridiagonal(dl[2:n-1], dd, du[1:n-2]), dl[1], du[end]
end

function pde_solve(prob::DiffusionProblem{R, T}, bc::ConstantCondition, 
        nx::Int, nt::Int) where {R<:Real, T<:Real}
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cdiff = dt / (2 * dx*dx) * prob.diffusivity
    println("Courant number is: $(abs(2 * cdiff))")
    if abs(2 * cdiff) > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    us = prob.initial.(xs)

    dv = one(cdiff) - 2*cdiff
    dd = fill(one(cdiff) + 2cdiff, nx-1)
    dl = fill(-cdiff, nx-2)
    matrixret = Tridiagonal(dl, dd, dl)
    bs = Vector{eltype(us)}(undef, nx-1)

    solution = Matrix{eltype(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us

    for i in 2:nt+1
        bs .= dv * us[2:nx] + cdiff * (us[1:nx-1] + us[3:end])
        bs[1] += cdiff * us[1]
        bs[end] += cdiff * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution
end

function pde_solve(prob::DiffusionProblem{F, T}, bc::ConstantCondition, 
        nx::Int, nt::Int) where {F<:Fct, T<:Real}
    
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cc = dt / (2 * dx*dx)
    println("Courant number is: $(2cc)")
    if 2cc > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    us = prob.initial.(xs)
    cdiff = cc * prob.diffusivity.(xs)

    matrix, edgel, edgeu = diffusionmatrix(cdiff)
    matrixret = I(nx-1) .- matrix
    matrixadv = I(nx-1) .+ matrix
    bs = Vector{eltype(us)}(undef, nx-1)

    solution = Matrix{eltype(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us

    for i in 2:nt+1
        bs .= matrixadv * us[2:nx]
        bs[1] += 2 * edgel * us[1]
        bs[end] += 2 * edgeu * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution
end

function totalmatrix(d::Vector{T}, c::Vector{T}) where T<:Real
    n = length(d) - 1
    dd = Vector{T}(undef, n-1)
    dl = Vector{T}(undef, n-1)
    du = Vector{T}(undef, n-1)
    for j in 1:n-1
        dd[j] = -2 * d[j+1] - c[j+2] + c[j] 
        dl[j] = d[j+1] - (d[j+2] - d[j]) / 4 + c[j+1]
        du[j] = c[j+1] + (c[j+2] - c[j]) / 4 - c[j+1]
    end
    return Tridiagonal(dl[2:n-1], dd, du[1:n-2]), dl[1], du[end]
end

function pde_solve(prob::ConstConvDiffProblem{R, T}, bc::ConstantCondition, nx::Int, 
        nt::Int) where {R<:Real, T<:Real}
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cdiff = dt / (2*dx*dx) * prob.diffusivity
    cconv = dt/ (4dx) * prob.velocity
    srate = dt * prob.sourcerate

    println("Courant number is: $(abs(2 * cdiff))")
    if abs(2 * cdiff) > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    us = prob.initial.(xs)

    dv = one(cdiff) - 2cdiff
    dd = fill(one(cdiff) + 2cdiff, nx-1)
    dp, dm = cdiff + cconv, cdiff - cconv
    dl = fill(-dp, nx - 2)
    du = fill(-dm, nx - 2)
    matrixret = Tridiagonal(dl, dd, du)
    bs = Vector{eltype(us)}(undef, nx-1)

    solution = Matrix{eltype(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us

    for i in 2:nt+1
        for j in 1:nx-1
            bs[j] = dv * us[j+1] + dp * us[j] + dm * us[j+2] + srate
        end
        bs[1] += dp * us[1]
        bs[end] += dm * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution
end

function pde_solve(prob::VariaConvDiffProblem{T}, bc::ConstantCondition, nx::Int, nt::Int) where T<:Real
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cd = dt / (2 * dx*dx)
    cc = dt / (4dx)
    println("Courant number is: $(2cd)")
    if 2cd > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2])

    xs = range(prob.xlim..., nx+1)
    us = prob.initial.(xs)

    cdiff = cd * prob.diffusivity.(xs)
    cconv = cc * prob.velocity.(xs)
    srate = dt * prob.sourcerate.(xs)

    matrix, edgel, edgeu = totalmatrix(cdiff, cconv)
    matrixret = I(nx-1) - matrix
    matrixadv = I(nx-1) + matrix
    bs = Vector{eltype(us)}(undef, nx-1)

    solution = Matrix{eltype(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us

    for i in 2:nt+1
        bs .= matrixadv * us[2:nx] + srate[2:nx]
        bs[1] += 2 * edgel * us[1]
        bs[end] += 2 * edgeu * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution
    
end

function pde_solve(prob::TimedDiffusionProblem{T}, bc::ConstantCondition, nx::Int, nt::Int) where T<:Real
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cc = dt / (2 * dx*dx)
    println("Courant number is: $(2cc)")
    if 2cc > 0.5
        @warn "Low numerical stability!"
    end

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2])

    xs = range(xlim..., nx+1)
    ts = range(T, tfin, nt+1)
    us = prob.initial.(xs)

    solution = Matrix{eltye(us)}(undef, nx+1, nt+1)
    solution[:, 1] .= us

    cdiff_new = cc .* prob.diffusivity.(xs, zero(dt))
    cdiff_old = Vector{eltype(cdiff_new)}(undef, nx+1)
    bs = Vector{eltype(us)}(undef, nx-1)

    for i in 2:nt+1
        cdiff_old .= cdiff_new
        cdiff_new .= prob.diffusivity.(xs, ts[i])

        matrixnew, edgel_new, edgeu_new = diffusionmatrix(cdiff_new)
        matrixold, edgel_old, edgeu_old = diffusionmatrix(cdiff_old)
        matrixret = I(nx-1) - matrixnew
        matrixadv = I(nx-1) + matrixold

        bs .= matrixadv * us[2:nx]
        bs[1] += (edgel_new + edgel_old) * us[1]
        bs[end] += (edgeu_new + edgeu_old) * us[end]

        us[2:nx] .= matrixret \ bs
        us[1] = left[2] - left[1] * us[2]
        us[end] = right[2] + right[1] * us[nx]

        solution[:, i] .= us
    end
    return solution
end