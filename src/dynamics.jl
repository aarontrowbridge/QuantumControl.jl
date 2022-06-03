module Dynamics

export MultiQubitSystem
export state_dim
export control_dim
export simulate

using LinearAlgebra
using StaticArrays
using ForwardDiff
using FiniteDiff

import RobotDynamics as RD

RD.@autodiff struct MultiQubitSystem <: RD.ContinuousDynamics
    nqubits::Int
    nstates::Int
    nqstates::Int
    isodim::Int
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
            nqubits, nstates, nqstates, isodim, H_drift, H_drive
        )
    end
end

RD.state_dim(model::MultiQubitSystem) = model.nstates + 3
RD.control_dim(::MultiQubitSystem) = 1

# for some reason, in-place dynamics is not working
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

function simulate(model, x0, U, t, dt, N)
    X = [x0]
    for i = 2:N
        push!(X, RD.discrete_dynamics(model, X[i-1], U[i-1], t[i-1], dt))
    end
    return X
end

end
