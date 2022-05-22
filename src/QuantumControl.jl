module QuantumControl


# export MultiQubitSystem
# # export state_dim
# # export control_dim
# # export dynamics!

# export py_sparse_matrix_to_julia
# export load_qutip_object
# export get_qutip_matrix


# # using Pkg
# # ENV["PYTHON"] = "/Users/aaron/miniconda3/envs/QOC/bin/python"
# # Pkg.build("PyCall")

# using PyCall

# using RobotDynamics
# using TrajectoryOptimization
# using LinearAlgebra
# using StaticArrays
# using SparseArrays
# using ForwardDiff
# using FiniteDiff

# import RobotDynamics as RD
# import TrajectoryOptimization as TO

using Reexport

include("dynamics.jl")
@reexport using .Dynamics

include("utils.jl")
@reexport using .Utils

end
