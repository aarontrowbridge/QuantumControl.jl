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
        if isodynamics
            H_drift_real = real(H_drift)
            H_drift_imag = imag(H_drift)
            H_drift_iso = SMatrix{isodim, isodim}(
                [H_drift_real -H_drift_imag;
                 H_drift_imag  H_drift_real]
            )

            H_drive_real = real(H_drive)
            H_drive_imag = imag(H_drive)
            H_drive_iso = SMatrix{isodim, isodim}(
                [H_drive_real -H_drive_imag;
                 H_drive_imag  H_drive_real]
            )

            return new{typeof(H_drift_iso)}(
                nqubits, nstates, nqstates, isodim, isodynamics, H_drift_iso, H_drive_iso
            )
        else
            return new{typeof(H_drift)}(
                nqubits, nstates, nqstates, isodim, isodynamics, H_drift, H_drive
            )
        end
    end
end

# state dim isomorphisms (C²)^(⨂ⁿ) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
RD.state_dim(model::MultiQubitSystem) = model.nstates + 3

# single control dimension
RD.control_dim(::MultiQubitSystem) = 1

# function RD.dynamics!(model::MultiQubitSystem, ẋ, x, u)
#     nstates = model.nstates

#     ctrl_inds = (nstates + 1):(nstates + 3)
#     ∫a, a, da = x[ctrl_inds]
#     dda = u[1]
#     xdot[ctrl_inds] .= [a, da, dda]

#     H = model.H_drift + a * model.H_drive

#     for i = 1:model.nqstates
#         ψ = getqstate(x, i, model.isodim)

#         ψ̇ = -im * H * ψ

#         setqstate!(ẋ, ψ̇, i, model.isodim)
#     end
#     return nothing
# end

function RD.dynamics(model::MultiQubitSystem, x, u)
    nstates = model.nstates

    ẋ = similar(x)

    ctrl_inds = (nstates + 1):(nstates + 3)
    ∫a, a, da = x[ctrl_inds]
    dda = u[1]
    ẋ[ctrl_inds] .= [a, da, dda]
    # ẋ[nstates + 1] = a
    # ẋ[nstates + 2] = da
    # ẋ[nstates + 3] = dda

    H = model.H_drift + a * model.H_drive

    schroedinger!(ẋ, x, H, model)

    return SVector{nstates+3}(ẋ)
end

const I2 = SMatrix{2,2}(I(2))
const im2 = SMatrix{2,2}([0 1; -1 0])

function schroedinger!(ẋ, x, H, model)
    Hreal = real(H)
    Himag = imag(H)
    for i = 1:model.nqstates
        ψiso = x[(1 + (i - 1) * model.isodim):(i * model.isodim)]
        ψ̇iso = (kron(I2, Himag) + kron(im2, Hreal)) * ψiso
        ẋ[(1 + (i - 1) * model.isodim):(i * model.isodim)] .= ψ̇iso
    end
    # for i = 1:model.nqstates
    #     ψ = getqstate(x, i, model.isodim)
    #     ψ̇ = -im * H * ψ
    #     setqstate!(ẋ, ψ̇, i, model.isodim)
    # end
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
