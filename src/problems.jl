module Problems

export SingleQubitProblem

using ..Dynamics
using ..QuantumLogic

using TrajectoryOptimization
using RobotDynamics
using Altro
using LinearAlgebra
using StaticArrays
using Random

const RD = RobotDynamics
const TO = TrajectoryOptimization

function SingleQubitProblem(
    H_drift,
    H_drive,
    gate::Symbol,
    ψ0::Union{Vector{T}, Vector{Vector{T}}} where T<:Number;
    kwargs...
)
    if typeof(ψ0) <: Vector{T} where T<:Number
        @assert size(H_drift)[2] == size(ψ0)[1] "Hamiltonian dimension does not match ket dimension"
        nqstates = 1
        ψf = apply(gate, ψ0)
        ψfs = [ψf]
        ψ̃f = ket_to_iso(ψf)
        ψ̃0 = ket_to_iso(ψ0)
    else
        @assert size(H_drift)[2] == size(ψ0[1])[1] "Hamiltonian dimension does not match ket dimension"
        nqstates = length(ψ0)
        ψfs = apply.(gate, ψ0)
        ψ̃fs = ket_to_iso.(ψfs)
        ψ̃f = foldr(vcat, ψ̃fs)
        ψ̃0s = ket_to_iso.(ψ0)
        ψ̃0 = foldr(vcat, ψ̃0s)
    end

    SingleQubitProblem(
        H_drift,
        H_drive,
        nqstates,
        ψ̃0,
        ψ̃f,
        ψfs;
        kwargs...
    )
end

function SingleQubitProblem(
    H_drift,                     # drift Hamiltonian
    H_drive,                     # drive Hamiltonian
    nqstates,                    # number of quantum states
    ψ̃0,                          # initial quantum state isomorphism(s)
    ψ̃f,                          # final quantum state isomorphism(s)
    ψfs;                         # vector final quantum states
    c0=zeros(3),                 # initial control input (∫a, a, ȧ) .= 0
    cf=zeros(3),                 # final control input (∫a, a, ȧ) .= 0
    dt=0.01,                     # time step
    N=1001,                      # number of knot points
    fidelity_cost=true,          # whether to use fidelity cost
    Qᵢᵢ=1.0,                     # diagonal cost for state
    Rᵢᵢ=0.1,                     # diagonal cost for control
    ctrl_cost=0.5,               # control cost
    state_cost_multiplier=100.0, # diagonal cost for final state
    ctrl_cost_multiplier=100.0,  # cost multiplier for state control variables
    dec_cost_multiplier=100.0,   # cost multiplier for decision variable
    U0=nothing,                  # initial control input
    U0_multiplier=1.0,           # initial random decision variable multiplier
    abound=true,                 # whether to bound the control input
    äbound=true,                 # whether to bound the decision variable
    äbound_val=2.0,              # bound value for decision variable
    agoal=true,                  # whether to constrain the control input to go to goal
    ψ̃goal=true,                  # whether to constrain the final state to be the goal state
)
    # final time
    tf = dt * (N - 1)

    # build model from Hamiltonians
    model = MultiQubitSystem(H_drift, H_drive, nqstates)
    n, m = RD.dims(model)

    # discretize model
    dmodel = RD.DiscretizedDynamics{RD.RK4}(model)

    # build initial and final state vectors
    x0 = SVector{n}([ψ̃0; c0])
    xf = SVector{n}([ψ̃f; cf])

    # random initial guess for control inputs
    if U0 === nothing
        U0 = [U0_multiplier * @SVector randn(m) for _ = 1:N-1]
    end

    if fidelity_cost
        Q = [@SVector fill(Qᵢᵢ, nqstates); @SVector fill(ctrl_cost, 3)]
        R = @SVector fill(Rᵢᵢ, 1)
        cost = MultiQubitSystemCost(ψfs; Q=Q, R=R)
        Qf = SVector{nqstates+3}([state_cost_multiplier * Q[1:end-3]; ctrl_cost_multiplier * Q[end-2:end]])
        Rf = dec_cost_multiplier * R
        Ri = Rf
        cost_term = MultiQubitSystemCost(ψfs; Q=Qf, R=Rf)
        obj = Objective(cost, cost_term, N)
        obj.cost[1] = MultiQubitSystemCost(ψfs; Q=Q, R=Ri)
    else
        # set up quadratic objective function
        Q = Diagonal(@SVector fill(Qᵢᵢ, n))
        R = Diagonal(@SVector fill(Rᵢᵢ, m))
        Qf = state_cost_multiplier * Q
        obj = LQRObjective(Q, R, Qf, xf, N; checks=false);

        # penalize control at first and last time step more
        costfun = LQRCost(Q, ctrl_cost_multiplier * R, xf)
        obj.cost[1] = costfun
        obj.cost[N-1] = costfun
    end

    # set up contraints
    cons = ConstraintList(n, m, N)

    # constraint: bound decision variable |ä| ≤ äbound_val

    ψ̃_max = ones(n-3)
    ψ̃_min = -ψ̃_max

    ∫a_max = Inf
    ∫a_min = -∫a_max

    a_max = 0.5
    a_min = -a_max

    ȧ_max = Inf
    ȧ_min = -ȧ_max

    ä_max = äbound ? äbound_val : Inf
    ä_min = -ä_max

    x_max = [ψ̃_max; ∫a_max; a_max; ȧ_max]
    x_min = [ψ̃_min; ∫a_min; a_min; ȧ_min]

    u_max = [ä_max]
    u_min = [ä_min]

    boundcon = BoundConstraint(n, m; x_min, x_max, u_min, u_max)
    add_constraint!(cons, boundcon, N)

    # if äbound
    #     äboundcon = NormConstraint(n, m, äbound_val, Altro.SecondOrderCone(), :control)
    #     add_constraint!(cons, äboundcon, N)
    # end

    # # constraint: bound control variable |aₖ| ≤ 0.5 GHz
    # if abound
    #     aboundcon = NormConstraint(n, m, 0.5, Altro.SecondOrderCone(), [model.nstates + 2])
    #     add_constraint!(cons, aboundcon, N)
    # end

    # constraint: set goal for ∫aₙ = aₙ = 0
    if agoal
        agoalcon = GoalConstraint(xf, model.nstates+1:model.nstates+2)
        add_constraint!(cons, agoalcon, N)
    end

    # constraint: set goal for |ψⁱₙ⟩ = |ψⁱfinal⟩
    if ψ̃goal
        # ψ̃goalcon = GoalConstraint(xf, 1:model.nstates)
        # add_constraint!(cons, ψ̃goalcon, N)

        ψgoalcon = QuantumGoalConstraint(ψfs)
        add_constraint!(cons, ψgoalcon, N)
    end

    return Problem(
        dmodel,
        obj,
        x0,
        tf;
        xf=xf,
        constraints=cons,
        U0=U0
    )
end

function Dynamics.simulate(prob::Problem, U)
    dmodel, _ = prob.model
    x0 = prob.x0
    N = prob.N
    t = gettimes(prob)
    dt = t[2] - t[1]
    X = simulate(dmodel, x0, U, t, dt, N)
    return X
end

end
