using QuantumControl
using Altro
using TrajectoryOptimization
using StaticArrays

println("\npackages are loaded!")

# getting Hamiltonians from saved qutip objects
qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

nqstates = 1

# X gate on single basis state: |0⟩ -> X |0⟩ = |1⟩
ψ̃00 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f0 = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

prob0 = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃00, ψ̃f0;
    N=1001,
)

solver0 = ALTROSolver(prob0; verbose=2)

# X gate on single basis state: |1⟩ -> X |1⟩ = |0⟩
ψ̃01 = @SVector [0.0, 1.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f1 = @SVector [1.0, 0.0, 0.0, 0.0] # final quantum state isomorphism

prob1 = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃01, ψ̃f1;
    N=1001,
)

solver1 = ALTROSolver(prob1; verbose=2)

println("\nbeginning ALTRO solves...")

println("\n|0⟩ -> X |0⟩")
solve!(solver0)

println("\n|1⟩ -> X |1⟩")
solve!(solver1)

println("\nplotting results...")

plot_wfn_ctrl_derivs(solver0, "plots/single_fluxonium_qubit/single_qstate_X_gate_on_0_basis_all.png")
plot_wfn_ctrl_derivs(solver1, "plots/single_fluxonium_qubit/single_qstate_X_gate_on_1_basis_all.png")

plot_wfn(solver0, "plots/single_fluxonium_qubit/single_qstate_X_gate_on_0_basis_wfn.png")
plot_wfn(solver1, "plots/single_fluxonium_qubit/single_qstate_X_gate_on_1_basis_wfn.png")
