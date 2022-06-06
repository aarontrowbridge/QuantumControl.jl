module QuantumControl

using Reexport

@reexport using RobotDynamics
@reexport using TrajectoryOptimization
@reexport using Altro
@reexport using StaticArrays
@reexport using LaTeXStrings

include("quantumlogic.jl")
@reexport using .QuantumLogic

include("dynamics.jl")
@reexport using .Dynamics

include("problems.jl")
@reexport using .Problems

include("qutiputils.jl")
@reexport using .QuTiPUtils

include("plottingutils.jl")
@reexport using .PlottingUtils

end
