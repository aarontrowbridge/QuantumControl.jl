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
    # qstates::SVector{N, SVector{M, T}} where {N, M, T}
    nqstates::Int
    H_drift::SMatrix{n} where {n}
    H_drive::SMatrix{m} where {m}

    function MultiQubitSystem(
        H_drift::SMatrix{n} where {n},
        H_drive::SMatrix{m} where {m},
        nqstates::Int
        )
        nqubits = Int(log2(size(H_drift)[1]))
        nstates = nqstates * 2^nqubits * 2
        return new{typeof(H_drift)}(
            nqubits, nstates, nqstates, H_drift, H_drive
        )
    end
end

# state dim isomorphisms (C²)^(⨂ⁿ) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
RD.state_dim(model::MultiQubitSystem) = model.nstates + 3

# single control dimension
RD.control_dim(::MultiQubitSystem) = 1

function RD.dynamics!(model::MultiQubitSystem, xdot, x, u)
    nstates = model.nstates

    ∫a, a, da = x[(nstates+1):(nstates+3)]
    dda = u[1]


    H = model.H_drift + a * model.H_drive

    for i = 1:model.nqstates
        ψ = getqstate(model, x, i)

        ψdot = -im * H * ψ

        setqstate!(model, xdot, ψdot, i)
    end

    xdot[(nstates+1):(nstates+3)] = [a, da, dda]

    return nothing
end

function getqstate(model::MultiQubitSystem, x, i)
    iso_dim = div(model.nstates, model.nqstates)
    ψreal = x[(1 + (i - 1) * iso_dim ):(i * div(iso_dim, 2))]
    ψimag = x[(div(iso_dim, 2) + 1):(i * iso_dim)]
    ψ = SA[(ψreal + im * ψimag)...]
    return ψ
end

function setqstate!(model::MultiQubitSystem, x, ψ, i)
    iso_dim = div(model.nstates, model.nqstates)
    for j in 1:div(iso_dim, 2)
        x[(i - 1) * iso_dim + j] = real(ψ[j])
        x[(i - 1) * iso_dim + div(iso_dim, 2) + j] = imag(ψ[j])
    end
end

end
