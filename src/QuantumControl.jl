module QuantumControl

using Reexport

@reexport using RobotDynamics
@reexport using TrajectoryOptimization
@reexport using Altro
@reexport using StaticArrays
@reexport using LaTeXStrings

include("dynamics.jl")
@reexport using .Dynamics

include("problems.jl")
@reexport using .Problems

include("utils.jl")
@reexport using .Utils

include("plottingutils.jl")
@reexport using .PlottingUtils

end
