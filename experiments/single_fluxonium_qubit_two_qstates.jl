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

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

println("\nsolving single quantum state prob...")
# X gate on single basis state |0⟩
nqstates = 1
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

initprob = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃0, ψ̃f;
    N=1001,
    U0_multiplier=0.1,
)

initsolver = ALTROSolver(initprob; verbose=2)

solve!(initsolver)

U0 = controls(initsolver)

println("\nsolving two quantum state problem (with U0 from init problem)...")

# X gate on two basis states |0⟩ and |1⟩
nqstates = 2
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f = @SVector [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] # final quantum state isomorphism

twoqstateprob = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃0, ψ̃f;
    N=1001,
    U0=U0,
)

solver = ALTROSolver(twoqstateprob; verbose=2)

solve!(solver)

println("\nplotting results...")

# plot_wfn_ctrl_derivs(solver, "plots/single_fluxonium_qubit.png")

# plot_wfn(solver, "plots/single_fluxonium_qubit_wfn.png")

plot_two_wfns(solver, "plots/single_fluxonium_qubit_two_wfns.png"; show_ctrl=true)
