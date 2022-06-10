module Dynamics

export MultiQubitSystem
export MultiQubitSystemCost
export QuantumGoalConstraint
export state_dim
export control_dim
export simulate
export ⊗
export Id2
export im2
export qubit_cost

using ..QuantumLogic

using TrajectoryOptimization
using RobotDynamics
using ForwardDiff
using FiniteDiff
using LinearAlgebra
using StaticArrays

const TO = TrajectoryOptimization
const RD = RobotDynamics

#
# multi qubit system model
#

RD.@autodiff struct MultiQubitSystem <: RD.ContinuousDynamics
    nqubits::Int
    isodim::Int
    nqstates::Int
    nstates::Int
    H_drift::SMatrix{n} where n
    H_drive::SMatrix{m} where m

    function MultiQubitSystem(
        H_drift::SMatrix{n} where n,
        H_drive::SMatrix{m} where m,
        nqstates::Int
    )
        nqubits = Int(log2(size(H_drift)[1]))
        isodim = 2 * 2^nqubits # (C²)^(⊗n) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
        nstates = nqstates * isodim
        return new{typeof(H_drift)}(
            nqubits, isodim, nqstates, nstates, H_drift, H_drive
        )
    end
end

RD.state_dim(model::MultiQubitSystem) = model.nstates + 3
RD.control_dim(::MultiQubitSystem) = 1

function RD.dynamics!(model::MultiQubitSystem, ẋ, x, u)
    ctrl_inds = (model.nstates + 1):(model.nstates + 3)
    ∫a, a, da = x[ctrl_inds]
    dda = u[1]
    ẋ[ctrl_inds] .= [a, da, dda]
    H = model.H_drift + a * model.H_drive
    schroedinger!(ẋ, x, H, model.nqstates, model.isodim)
end

function RD.dynamics(model::MultiQubitSystem, x, u)
    ẋ = similar(x)
    ctrl_inds = (model.nstates + 1):(model.nstates + 3)
    ∫a, a, ȧ = x[ctrl_inds]
    ä = u[1]
    ẋ[ctrl_inds] .= [a, ȧ, ä]
    H = model.H_drift + a * model.H_drive
    schroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    return ẋ
end

const Id2 = SMatrix{2,2}(I(2))
const im2 = SMatrix{2,2}([0 -1; 1 0])

⊗(A, B) = kron(A, B)

function schroedinger!(ẋ, x, H, nqstates, isodim)
    Hreal = real(H)
    Himag = imag(H)
    opr = Id2 ⊗ Himag - im2 ⊗ Hreal # this is isomorphism of -iH
    for i = 1:nqstates
        ψ̃ = x[(1 + (i - 1) * isodim):(i * isodim)]
        ψ̃dot = opr * ψ̃
        ẋ[(1 + (i - 1) * isodim):(i * isodim)] .= ψ̃dot
    end
end

function simulate(dmodel, x0, U, t, dt, N)
    X = [x0]
    for i = 2:N
        push!(X, RD.discrete_dynamics(dmodel, X[i-1], U[i-1], t[i-1], dt))
    end
    return X
end

#
# quantum cost function
#

struct MultiQubitSystemCost <: TO.CostFunction
    Q::SVector{N} where N
    R::SVector{M} where M
    ψfs::Vector{SVector{D, T}} where {D,T}
    nqstates::Int
    isodim::Int
    function MultiQubitSystemCost(
        ψf;
        Q=[ones(isa(ψf, Vector{V} where V<:AbstractArray) ? length(ψf) : 1); ones(3)],
        R=[0.1],
    )
        if isa(ψf, Vector{A} where A<:AbstractArray)
            nqstates = length(ψf)
            isodim = 2 * length(ψf[1])
            ψfs = [ket(ψ) for ψ in ψf]
        else
            nqstates = 1
            isodim = 2 * length(ψf)
            ψfs = [ket(ψf)]
        end
        return new(Q, R, ψfs, nqstates, isodim)
    end
end

RD.state_dim(cost::MultiQubitSystemCost) = cost.isodim * cost.nqstates + 3
RD.control_dim(::MultiQubitSystemCost) = 1

(cost::MultiQubitSystemCost)(xu) = RD.evaluate(cost, xu[1:end-1], xu[end:end])

function RD.evaluate(cost::MultiQubitSystemCost, x, u)
    ψ̃s = [x[(1 + (i-1)*cost.isodim):i*cost.isodim] for i in 1:cost.nqstates]
    ψs = iso_to_ket.(ψ̃s)
    J = 0.0
    for (i, (ψ, ψf)) in enumerate(zip(ψs, cost.ψfs))
        J += cost.Q[i] * qubit_cost(ψ, ψf)
        # @info qubit_cost(ψ, ψf)
    end
    J += dot(cost.Q[end-2:end], x[end-2:end].^2)
    # @info dot(cost.Q[end-2:end], x[end-2:end].^2)
    if !isempty(u)
        J += cost.R[1] * u[1]^2
        # @info cost.R[1] * u[1]^2
    end
    if J == NaN exit() end
    @info "J = " J
    return J
end

qubit_cost(ψ, ψf) = min(abs(1 - ψ'ψf), abs(1 + ψ'ψf))

function RD.gradient!(cost::MultiQubitSystemCost, grad, x, u)
    ForwardDiff.gradient!(grad, cost, [x; !isempty(u) ? u : 0])
end

function RD.hessian!(cost::MultiQubitSystemCost, hess, x, u)
    ForwardDiff.hessian!(hess, cost, [x; !isempty(u) ? u : 0])
end

#
# quantum goal constraint
#

struct QuantumGoalConstraint <: TO.StateConstraint
    ψfs::Vector{SVector{N,T}} where {N,T<:Number}
    nqstates::Int
    isodim::Int

    function QuantumGoalConstraint(
        ψfs::Vector{SVector{N,T}} ,
    ) where {N,T<:Number}
        nqstates = length(ψfs)
        isodim = 2 * N
        return new(ψfs, nqstates, isodim)
    end
end

TO.sense(::QuantumGoalConstraint) = Equality()

RD.state_dim(con::QuantumGoalConstraint) = con.isodim * con.nqstates + 3
RD.control_dim(::QuantumGoalConstraint) = 1
RD.output_dim(::QuantumGoalConstraint) = 1

(con::QuantumGoalConstraint)(x) = RD.evaluate(con, x)

RD.evaluate(con::QuantumGoalConstraint, z::RD.AbstractKnotPoint) = con(RD.state(z))
RD.evaluate!(con::QuantumGoalConstraint, c, x) = c .= con(x)

function RD.evaluate(con::QuantumGoalConstraint, x)
    ψ̃s = [x[(1 + (i-1)*con.isodim):i*con.isodim] for i = 1:con.nqstates]
    ψs = iso_to_ket.(ψ̃s)
    c = zeros(typeof(x[1]), 1)
    for (ψ, ψf) in zip(ψs, con.ψfs)
        c .+= qubit_cost(ψ, ψf)
    end
    return c
end

RD.jacobian!(con::QuantumGoalConstraint, jac, x) = ForwardDiff.jacobian!(jac, con, x)
RD.jacobian!(::RD.StateOnly, con::QuantumGoalConstraint, jac, y, z::RD.AbstractKnotPoint) =
    RD.jacobian!(con, jac, RD.state(z))

end
