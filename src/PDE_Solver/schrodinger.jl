include("utils.jl")

mutable struct StaticSchrodingerProblem{T <: Real} <: AbstractProblem
    mass::Number
    potential::Fct
    inital_wavefunc::Fct
    xlim::Tuple{T, T}
    tfin::T
end

mutable struct DynamicSchrodingerProblem{T <: Real} <: AbstractProblem
    mass::T
    potential::Fct
    inital_wavefunc::Fct
    xlim::Tuple{T, T}
    tfin::T
end

mutable struct EMSchrodingerProblem{T <: Real} <: AbstractProblem
    mass::T
    charge::T
    ext_potential::Fct
    mag_potential::Fct
    ele_potential::Fct
    inital_wavefunc::Fct
    xlim::Tuple{T, T}
    tfin::T
end


function diagval(c::Number, v::Real, dt::Real)
    return  2*c + complex(0, dt * v / 2)
end

function diagval(c0::Number, c1::Real, c2::Number, v::Real, a::Vector{<:Real}, qp::Real, dt)
    return 2*c0 + c2*a[2]*a[2] + complex(0 - c1*(a[3] - a[1])/2, qp*dt/2 + v*dt/2)
end

function pde_solve(prob::StaticSchrodingerProblem, bc::ConstantCondition, nx::Int, nt::Int)
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cc = im * dt / (4 * prob.mass * dx*dx)
    println("Courant number is: $(abs(cc))")
    if abs(cc) > 0.5
        @warn "Low numerical stability!"
    end 

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    psi = complex.(prob.inital_wavefunc.(xs))

    diagp = zeros(typeof(cc), nx - 1)
    diagl  = fill(cc, nx - 2)
    for (j, x) in enumerate(xs[2:nx])
        v = prob.potential(x)
        diagp[j] -= diagval(cc, v, dt)
    end
    
    matrix = Tridiagonal(diagl, diagp, diagl)
    matrixret = I(nx-1) .- matrix
    matrixadv = I(nx-1) .+ matrix
    bs = Vector{eltype(psi)}(undef, nx-1)

    solution = Matrix{eltype(psi)}(undef, nx+1, nt+1)
    solution[:, 1] .= psi

    for i in 2:(nt + 1)
        bs .= matrixadv * psi[2:nx]
        bs[1] += 2 * cc * psi[1]
        bs[end] += 2 * cc * psi[end]

        psi[2:nx] .= matrixret \ bs
        psi[1] = left[2] - left[1] * psi[2]
        psi[end] = right[2] + right[1] * psi[nx]

        solution[:, i] .= psi
    end
    return solution
end

function pde_solve(prob::EMSchrodingerProblem, bc::ConstantCondition, nx::Int, nt::Int)
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cc = im * dt / (4 * prob.mass * dx*dx)
    q = prob.charge
    c1 = q  / prob.mass * dt / dx / 4
    c2 = im * (q/ prob.mass)^2 * dt / 2
    println("Courant number is: $(abs(cc))")
    if abs(cc) > 0.5
        @warn "Low numerical stability!"
    end 

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)

    avec = prob.avec.(xs)
    psi = complex.(prob.inital_wavefunc.(xs))

    diagp = Vector{typeof(cc)}(undef, nx-1)
    diagl = fill(cc, nx - 2)
    diagu = fill(cc, nx - 2)
    for (j, x) in enumerate(xs[2:nx])
        aloc = avec[j:j+2]
        vext = prob.ext_potential(x)
        phie = prob.ele_potential(x)

        diagp[j] = - diagval(cc, c1, c2, vext, aloc, q*phie, dt)
        if j < nx-1
            diagl[j] -= c1*aloc[2]
            diagu[j] += c1*aloc[2]
        end
    end

    matrix = Tridiagonal(diagl, diagp, diagu)
    matrixret = I(nx-1) .- matrix
    matrixadv = I(nx-1) .+ matrix
    bs = Vector{eltype(psi)}(undef, nx-1)

    solution = Matrix{eltype(psi)}(undef, nx+1, nt+1)
    solution[:, 1] .= psi

    for i in 2:nt+1
        bs .= matrixadv * psi[2:nx]
        bs[1] += 2 * (cc - c1*avec[2]) * psi[1]
        bs[end] += 2 * (cc + c1*avec[nx]) * psi[end]

        psi[2:nx] .= matrixret \ bs
        psi[1] = left[2] - left[1] * psi[2]
        psi[end] = right[2] + right[1] * psi[nx]

        solution[:, i] .= psi
    end
    return solution
end

function pde_solve(prob::DynamicSchrodingerProblem, bc::ConstantCondition, nx::Int, nt::Int)
    dx = (prob.xlim[2] - prob.xlim[1]) / nx
    dt = prob.tfin / nt
    cc = im * dt / (4 * prob.mass * dx*dx)
    println("Courant number is: $(abs(cc))")
    if abs(cc) > 0.5
        @warn "Low numerical stability!"
    end 

    left = (bc.bcl[2], bc.bcl[3] * dx) ./ (bc.bcl[1] * dx - bc.bcl[2])
    right = (bc.bch[2], bc.bch[3] * dx) ./ (bc.bch[1] * dx + bc.bch[2]) 

    xs = range(prob.xlim..., nx+1)
    ts = range(zero(dt), prob.tfin, nt+1)
    psi = complex.(prob.inital_wavefunc.(xs))

    diagpp = ones(typeof(cc), nx - 1)
    diagl = fill(-cc, nx - 2)
    matrix = Tridiagonal(diagl, diagpp, diagl)
    bs = Vector{eltype(psi)}(undef, nx-1)

    solution = Matrix{eltype(psi)}(undef, nx+1, nt+1)
    solution[:, 1] .= psi

    potential_new = prob.potential.(xs[2:nx], zero(dt))
    potential_old = Vector{eltype(potential_new)}(undef, nx-1)
    for i in 2:(nt + 1)
        potential_old .= potential_new
        potential_new .= prob.potential.(xs[2:nx], ts[i])

        for j in 1:nx-1
            matrix[j, j] = one(cc) + diagval(cc, potential_new[j], dt)
            diagpm = one(cc) - diagval(cc, potential_old[j-1], dt)
            bs[j] = diagpm * psi[j+1] + cc * (psi[j] + psi[j+2])
        end
        bs[1] += cc * psi[1]
        bs[end] += cc * psi[end] 

        psi[2:nx] .= matrix \ bs
        psi[1] = left[2] - left[1] * psi[2]
        psi[end] = right[2] + right[1] * psi[nx]

        solution[:, i] .= psi
    end
    return solution
end

