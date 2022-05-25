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
    nqstates::Int
    iso_dim::Int
    H_drift::SMatrix{n} where n
    H_drive::SMatrix{m} where m

    function MultiQubitSystem(
        H_drift::SMatrix{n} where n,
        H_drive::SMatrix{m} where m,
        nqstates::Int
        )
        nqubits = Int(log2(size(H_drift)[1]))
        iso_dim = 2^nqubits * 2
        nstates = nqstates * iso_dim
        return new{typeof(H_drift)}(
            nqubits, nstates, nqstates, iso_dim, H_drift, H_drive
        )
    end
end

# state dim isomorphisms (C²)^(⨂ⁿ) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
RD.state_dim(model::MultiQubitSystem) = model.nstates + 3

# single control dimension
RD.control_dim(::MultiQubitSystem) = 1

function RD.dynamics!(model::MultiQubitSystem, ẋ, x, u)
    nstates = model.nstates

    ∫a, a, da = x[(nstates+1):(nstates+3)]
    dda = u[1]


    H = model.H_drift + a * model.H_drive

    for i = 1:model.nqstates
        ψ = getqstate(x, i, model.iso_dim)

        ψ̇ = -im * H * ψ

        setqstate!(ẋ, ψ̇, i, model.iso_dim)
    end

    xdot[(nstates+1):(nstates+3)] = [a, da, dda]

    return nothing
end

function RD.dynamics(model::MultiQubitSystem, x, u)
    nstates = model.nstates

    ẋ = @MVector zeros(nstates + 3)

    ∫a, a, da = x[(nstates+1):(nstates+3)]
    dda = u[1]


    H = model.H_drift + a * model.H_drive

    for i = 1:model.nqstates
        ψ = getqstate(x, i, model.iso_dim)

        ψ̇ = -im * H * ψ

        setqstate!(ẋ, ψ̇, i, model.iso_dim)
    end

    ẋ[(nstates+1):(nstates+3)] .= [a, da, dda]

    return SVector{nstates+3}(ẋ)
end


function getqstate(x, i, iso_dim)
    ψreal = x[(1 + (i - 1) * iso_dim ):(i * div(iso_dim, 2))]
    ψimag = x[(div(iso_dim, 2) + 1):(i * iso_dim)]
    ψ = SA[(ψreal + im * ψimag)...]
    return ψ
end

function setqstate!(x, ψ, i, iso_dim)
    real_inds = (1 + (i - 1) * iso_dim):(div(iso_dim, 2))
    imag_inds = real_inds .+ div(iso_dim, 2)
    x[real_inds] .= real(ψ)
    x[imag_inds] .= imag(ψ)
end

end
