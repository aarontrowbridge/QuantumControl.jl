module QuantumLogic

export apply
export ket_to_iso
export iso_to_ket

using StaticArrays

const GATES = Dict(
    :X => [0 1;
           1 0],

    :Y => [0 -im;
           im 0],

    :Z => [1 0;
           0 -1],

    :H => [1 1;
           1 -1]/√2,

    :CX => [1 0 0 0;
            0 1 0 0;
            0 0 0 1;
            0 0 1 0],
)

function apply(gate::Symbol, ψ)
    @assert gate in keys(GATES) "gate not found"
    U = GATES[gate]
    @assert size(U)[2] == size(ψ)[1] "gate size does not match ket dim"
    return  U * ψ
end

function ket_to_iso(ψ)
    isodim = 2 * size(ψ)[1]
    ψreal = real(ψ)
    ψimag = imag(ψ)
    ψ̃ = [ψreal; ψimag]
    return SVector{isodim, Float64}(ψ̃)
end

function iso_to_ket(ψ̃)
    ketdim = div(size(ψ̃)[1], 2)
    ψreal = ψ̃[1:ketdim]
    ψimag = ψ̃[ketdim+1:2*ketdim]
    ψ = ψreal + im * ψimag
    return ψ
end


end
