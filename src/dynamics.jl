module Dynamics

export MultiQubitSystem
export state_dim
export control_dim
export dynamics!


using RobotDynamics
using TrajectoryOptimization
using LinearAlgebra
using StaticArrays
using ForwardDiff
using FiniteDiff

import RobotDynamics as RD
import TrajectoryOptimization as TO

RD.@autodiff struct MultiQubitSystem <: RD.ContinuousDynamics
    nqubits::Int
    nstates::Int
    H_drift::SMatrix{n} where {n}
    H_drive::SMatrix{m} where {m}

    function MultiQubitSystem(
        H_drift::SMatrix{n} where {n},
        H_drive::SMatrix{m} where {m}
        )
        nqubits = Int(log2(size(H_drift)[1]))
        nstates = 2^nqubits * 2
        return new{typeof(H_drift)}(
            nqubits, nstates, H_drift, H_drive
        )
    end
end

# state dim isomorphisms (C²)^(⨂ⁿ) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
RD.state_dim(model::MultiQubitSystem) = model.nstates + 3

# single control dimension
RD.control_dim(::MultiQubitSystem) = 1

function RD.dynamics!(model::MultiQubitSystem, xdot, x, u)
    nstates = nstates(model)

    ∫a, a, da = x[(nstates+1):(nstates+3)]
    dda = u[1]


    H = model.H_drift + a * model.H_drive

    ψ = getqstate(model, x)

    ψdot = -im * H * ψ

    setqstate!(model, xdot, ψdot)

    xdot[(nstates+1):(nstates+3)] = [a, da, dda]

    return nothing
end

function getqstate(model::MultiQubitSystem, x)
    nstates = model.nstates
    ψreal = x[1:div(nstates, 2)]
    ψimag = x[div(nstates, 2)+1:nstates]
    ψ = SA[(ψreal + im * ψimag)...]
    return ψ
end

function setqstate!(model::MultiQubitSystem, x, ψ)
    nstates = model.nstates
    for i in 1:div(nstates, 2)
        x[i] = real(ψ[i])
        x[nstates + i] = imag(ψ[i])
    end
end

end
