module Dynamics

export MultiQubitSystem
export state_dim
export control_dim

using LinearAlgebra
using StaticArrays
using ForwardDiff
using FiniteDiff

import RobotDynamics as RD

RD.@autodiff struct MultiQubitSystem <: RD.ContinuousDynamics
    nqubits::Int
    nstates::Int
    nqstates::Int
    ctrl_order::Int
    isodim::Int
    isodynamics::Bool
    H_drift::SMatrix{n} where n
    H_drive::SMatrix{m} where m

    function MultiQubitSystem(
        H_drift::SMatrix{n} where n,
        H_drive::SMatrix{m} where m,
        nqstates::Int,
        ctrl_order::Int;
        isodynamics=true
    )
        nqubits = Int(log2(size(H_drift)[1]))
        isodim = 2 * 2^nqubits # (C²)^(⊗n) ≅ C^(2ⁿ) ≅ R^(2⋅2ⁿ)
        nstates = nqstates * isodim
        return new{typeof(H_drift)}(
            nqubits, nstates, nqstates, ctrl_order, isodim, isodynamics, H_drift, H_drive
        )
    end
end

RD.state_dim(model::MultiQubitSystem) = model.nstates + model.ctrl_order + 1
RD.control_dim(::MultiQubitSystem) = 1

# for some reason, in-place dynamics is not working
function RD.dynamics!(model::MultiQubitSystem, ẋ, x, u)
    ctrl_inds = (model.nstates + 1):(model.nstates + model.ctrl_order + 1)
    ∫a, a, da = x[ctrl_inds]
    dda = u[1]
    ẋ[ctrl_inds] .= [a, da, dda]
    H = model.H_drift + a * model.H_drive
    if model.isodynamics
        isoschroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    else
        schroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    end
end

function RD.dynamics(model::MultiQubitSystem, x, u)
    ẋ = similar(x)
    ctrl_inds = (model.nstates + 1):(model.nstates + model.ctrl_order + 1)
    ∫a, a, ȧ = x[ctrl_inds]
    ä = u[1]
    ẋ[ctrl_inds] .= [a, ȧ, ä]
    H = model.H_drift + a * model.H_drive
    if model.isodynamics
        isoschroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    else
        schroedinger!(ẋ, x, H, model.nqstates, model.isodim)
    end
    return ẋ
end

const Id2 = SMatrix{2,2}(I(2))
const im2 = SMatrix{2,2}([0 -1; 1 0])

⊗(A, B) = kron(A, B)

function isoschroedinger!(ẋ, x, H, nqstates, isodim)
    Hreal = real(H)
    Himag = imag(H)
    opr = Id2 ⊗ Himag - im2 ⊗ Hreal # this is isomorphism of -iH
    for i = 1:nqstates
        ψ̃ = x[(1 + (i - 1) * isodim):(i * isodim)]
        ψ̃dot = opr * ψ̃
        ẋ[(1 + (i - 1) * isodim):(i * isodim)] .= ψ̃dot
    end
end

#
# below, direct schroedinger does not work with autodiff
# need to implment RD.jacobian!
#

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
