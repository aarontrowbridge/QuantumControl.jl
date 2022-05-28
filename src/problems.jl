module Problems

export SingleQubitProblem

using ..Dynamics

using TrajectoryOptimization
using Altro
using LinearAlgebra
using StaticArrays
using Random

import RobotDynamics as RD

function SingleQubitProblem(
    H_drift,                    # drift Hamiltonian
    H_drive,                    # drive Hamiltonian
    nqstates,                   # number of quantum states
    ψ̃0,                         # initial quantum state isomorphism
    ψ̃f;                         # final quantum state isomorphism
    order=2,                    # derivative order of decision variable (dda)
    c0=zeros(order+1),          # initial control input (∫a, a, ȧ) .= 0
    cf=zeros(order+1),          # final control input (∫a, a, ȧ) .= 0
    dt=0.01,                    # time step
    N=101,                      # number of knot points
    Qᵢᵢ=1.0,                    # diagonal cost for state
    Rᵢᵢ=0.1,                    # diagonal cost for control
    Qfᵢᵢ=100.0,                 # diagonal cost for final state
    ctrl_cost_multiplier=100.0, # cost multiplier for control
    U0_multiplier=1.0,          # initial random decision variable multiplier
)
    # final time
    tf = dt * (N - 1)

    # build model from Hamiltonians
    model = MultiQubitSystem(H_drift, H_drive, nqstates, order)
    n, m = RD.dims(model)

    # discretize model
    dmodel = RD.DiscretizedDynamics{RD.RK4}(model)

    # build initial and final state vectors
    x0 = SVector{n}([ψ̃0; c0])
    xf = SVector{n}([ψ̃f; cf])

    # random initial guess for control inputs
    U0 = [U0_multiplier * @SVector randn(m) for _ = 1:N-1]

    # set up quadratic objective function
    Q = Diagonal(@SVector fill(Qᵢᵢ, n))
    R = Diagonal(@SVector fill(Rᵢᵢ, m))
    Qf = Diagonal(@SVector fill(Qfᵢᵢ, n))
    obj = LQRObjective(Q, R, Qf, xf, N; checks=false);

    # penalize control at first and last time step more
    costfun = LQRCost(Q, ctrl_cost_multiplier * R, xf)
    obj.cost[1] = costfun
    obj.cost[N-1] = costfun

    # set up contraints
    cons = ConstraintList(n, m, N)

    # constraint: maintain normalized quantum statevectors
    normcon = NormConstraint(n, m, 1.0*model.nqstates, Equality(), 1:model.nstates)
    add_constraint!(cons, normcon, N)

    # constraint: bound control variable |aₖ| ≤ 0.5 GHz
    aboundcon = NormConstraint(n, m, 0.5e9, Inequality(), [model.nstates + 2])
    add_constraint!(cons, aboundcon, N)

    # # constraint: set goal for ∫a₁ = a₁ = ȧ₁ = 0
    # a₁goalcon = GoalConstraint(x0, model.nstates+1:model.nstates+order+1)
    # add_constraint!(cons, a₁goalcon, N)

    # constraint: set goal for ∫aₙ = aₙ = 0
    aₙgoalcon = GoalConstraint(xf, model.nstates+1:model.nstates+order)
    add_constraint!(cons, aₙgoalcon, N)

    # constraint: set goal for |ψⁱₙ⟩ = |ψⁱfinal⟩
    ψgoalcon = GoalConstraint(xf, 1:model.nstates)
    add_constraint!(cons, ψgoalcon, N)

    return Problem(dmodel, obj, x0, tf; xf=xf, constraints=cons, U0=U0)
end






end
