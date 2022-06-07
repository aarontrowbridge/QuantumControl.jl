module Dynamics

export MultiQubitSystem
export MultiQubitSystemCost
export state_dim
export control_dim
export simulate

using ..QuantumLogic

using TrajectoryOptimization
using RobotDynamics
using ForwardDiff
using FiniteDiff
using LinearAlgebra
using StaticArrays

const TO = TrajectoryOptimization
const RD = RobotDynamics

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

struct MultiQubitSystemCost <: TO.CostFunction
    Q::SVector{N} where N
    R::SVector{M} where M
    ψ̃f::SVector{L} where L
    nqstates::Int
    isodim::Int
end

(cost::MultiQubitSystemCost)(x, u) = RD.evaluate(cost, x, u)
(cost::MultiQubitSystemCost)(xu) = cost(xu[1:end-1], xu[end:end])

RD.state_dim(cost::MultiQubitSystemCost) = cost.isodim * cost.nqstates + 3
RD.control_dim(::MultiQubitSystemCost) = 1

function RD.evaluate(cost::MultiQubitSystemCost, x, u)
    ψ̃s = [x[(1 + (i-1)*cost.isodim):i*cost.isodim] for i in 1:cost.nqstates]
    ψs = iso_to_ket.(ψ̃s)
    ψ̃fs = [cost.ψ̃f[(1 + (i-1)*cost.isodim):i*cost.isodim] for i in 1:cost.nqstates]
    ψfs = iso_to_ket.(ψ̃fs)
    J = 0.0
    for (i, (ψ, ψf)) in enumerate(zip(ψs, ψfs))
        J += cost.Q[i] * (1 - abs2(ψ'ψf))^2
    end
    J += dot(cost.Q[end-2:end], x[end-2:end].^2)
    if !isempty(u)
        J += cost.R[1] * u[1]^2
    end
    return J
end

function RD.gradient!(cost::MultiQubitSystemCost, grad, x, u)
    if !isempty(u)
        ForwardDiff.gradient!(grad, cost, [x; u])
    else
        ForwardDiff.gradient!(grad, cost, [x; 0.0])
    end
end

function RD.hessian!(cost::MultiQubitSystemCost, hess, x, u)
    if !isempty(u)
        ForwardDiff.hessian!(hess, cost, [x; u])
    else
        ForwardDiff.hessian!(hess, cost, [x; 0.0])
    end
end

end
