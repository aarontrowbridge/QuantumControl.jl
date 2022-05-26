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
    isodim::Int
    isodynamics::Bool
    H_drift::SMatrix{n} where n
    H_drive::SMatrix{m} where m

    function MultiQubitSystem(
        H_drift::SMatrix{n} where n,
        H_drive::SMatrix{m} where m,
        nqstates::Int;
        isodynamics=true
    )
        nqubits = Int(log2(size(H_drift)[1]))
        isodim = 2^nqubits * 2
        nstates = nqstates * isodim
        return new{typeof(H_drift)}(
            nqubits, nstates, nqstates, isodim, isodynamics, H_drift, H_drive
        )
    end
end

# state dim isomorphisms (C²)^(⨂ⁿ) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
RD.state_dim(model::MultiQubitSystem) = model.nstates + 3

# single control dimension
RD.control_dim(::MultiQubitSystem) = 1

function RD.dynamics!(model::MultiQubitSystem, ẋ, x, u)
    nstates = model.nstates

    ctrl_inds = (nstates + 1):(nstates + 3)
    ∫a, a, da = x[ctrl_inds]
    dda = u[1]

    xdot[ctrl_inds] .= [a, da, dda]

    H = model.H_drift + a * model.H_drive

    if model.isodynamics
        isoschroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    else
        schroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    end
end

function RD.dynamics(model::MultiQubitSystem, x, u)
    ẋ = similar(x)

    nstates = model.nstates

    ctrl_inds = (nstates + 1):(nstates + 3)
    ∫a, a, da = x[ctrl_inds]
    dda = u[1]

    ẋ[ctrl_inds] .= [a, da, dda]

    H = model.H_drift + a * model.H_drive

    if model.isodynamics
        isoschroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    else
        schroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    end

    return SVector{nstates+3}(ẋ)
end

const Id2 = SMatrix{2,2}(I(2))
const im2 = SMatrix{2,2}([0 1; -1 0])

function isoschroedinger!(ẋ, x, H, nqstates, isodim)
    Hreal = real(H)
    Himag = imag(H)
    for i = 1:nqstates
        ψiso = x[(1 + (i - 1) * isodim):(i * isodim)]
        ψ̇iso = (kron(Id2, Himag) + kron(im2, Hreal)) * ψiso
        ẋ[(1 + (i - 1) * isodim):(i * isodim)] .= ψ̇iso
    end
end

function schroedinger!(ẋ, x, H, nqstates, isodim)
    for i = 1:nqstates
        ψ = getqstate(x, i, isodim)
        ψ̇ = -im * H * ψ
        setqstate!(ẋ, ψ̇, i, isodim)
    end
end

function getqstate(x, i, isodim)
    ψreal = x[(1 + (i - 1) * isodim ):(i * div(isodim, 2))]
    ψimag = x[(div(isodim, 2) + 1):(i * isodim)]
    ψ = SA[(ψreal + im * ψimag)...]
    return ψ
end

function setqstate!(x, ψ, i, isodim)
    real_inds = (1 + (i - 1) * isodim):(div(isodim, 2))
    imag_inds = real_inds .+ div(isodim, 2)
    x[real_inds] .= real(ψ)
    x[imag_inds] .= imag(ψ)
end

end
