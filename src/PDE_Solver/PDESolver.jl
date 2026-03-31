module PDESolver

using MKL, LinearAlgebra

include("utils.jl")
include("convection_diffusion.jl")
include("schrodinger.jl")
export ConvectionProblem, DiffusionProblem, ConstConvDiffProblem, VariaConvDiffProblem, TimedDiffusionProblem
export StaticSchrodingerProblem, DynamicSchrodingerProblem, EMSchrodingerProblem, pde_solve
end # module PDESolver