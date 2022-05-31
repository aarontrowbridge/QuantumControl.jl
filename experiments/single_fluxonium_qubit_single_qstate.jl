using QuantumControl
using Altro
using TrajectoryOptimization
using RobotDynamics
using LinearAlgebra
using StaticArrays
using Random
using CairoMakie

import TrajectoryOptimization as TO
import RobotDynamics as RD

println("\npackages are loaded!")

# getting Hamiltonians from saved qutip objects
qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

# X gate on single basis state: |0⟩ -> X |0⟩ = |1⟩
nqstates = 1
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

prob = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃0, ψ̃f;
    N=1001,
)

solver = ALTROSolver(prob; verbose=2)

println("\nbeginning ALTRO solve...")

solve!(solver)

println("\nplotting results...")

plot_wfn_ctrl_derivs(solver, "plots/single_fluxonium_qubit.png")

# plot_wfn(solver, "plots/single_fluxonium_qubit_wfn.png")
