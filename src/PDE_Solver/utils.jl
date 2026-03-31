using MKL
using LinearAlgebra

abstract type AbstractProblem end
abstract type AbstractCondition end

const Tuple3{T} = Tuple{T, T, T}
const Tuple3F = Tuple{Function, Function, Function}
const Fct = Function

function validcheck(x::Tuple3{T})::Bool where T <: Number 
    if iszero(x[1]) && iszero(x[2])
        return false
    else
        return true
    end
end

struct ConstantCondition{T <: Number} <: AbstractCondition
    bcl::Tuple3{T}
    bch::Tuple3{T}

    ConstantCondition{T}(bcl::Tuple, bch::Tuple) where T<:Number = new{T}(bcl, bch)
    ConstantCondition(bcl::Tuple3{T}, bch::Tuple3{T}) where T<:Number = 
        validcheck(bcl) && validcheck(bch) ? new{T}(bcl, bch) : error("Invalid boundary condition!")
end

struct VariantCondition <: AbstractCondition
    bcl::Tuple3F
    bch::Tuple3F
end

DirichletCondition(bcl::T, bch::T) where T <: Number = 
    ConstantCondition{T}((one(T), zero(T), bcl), (one(T), zero(T), bch))
DirichletCondition(bcl::Fct, bch::Fct) = 
    VariantCondition((one, zero, bcl), (one, zero, bch))

NeumannCondition(bcl::T, bch::T) where T <: Number = 
    ConstantCondition{T}((zero(T), one(T), bcl), (zero(T), one(T), bch))
NeumannCondition(bcl::Fct, bch::Fct) = 
    VariantCondition((zero, one, bcl), (zero, one, bch))