using QuantumControl

println("\npackages are loaded!")

# getting Hamiltonians from saved qutip objects
qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

N = 1001
nqstates = 1

# qubit basis states |0⟩ and |1⟩
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃1 = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

# X gate on single basis state: |0⟩ -> X |0⟩ = |1⟩
prob0 = SingleQubitProblem(
    H_drift,
    H_drive,
    nqstates,
    ψ̃0, # initial state
    ψ̃1; # final state
    N=N,
)

solver0 = ALTROSolver(prob0; verbose=2)

# X gate on single basis state: |1⟩ -> X |1⟩ = |0⟩
prob1 = SingleQubitProblem(
    H_drift,
    H_drive,
    nqstates,
    ψ̃1, # initial state
    ψ̃0; # final state
    N=N,
)

solver1 = ALTROSolver(prob1; verbose=2)

println("\nbeginning ALTRO solves...")

println("\n|0⟩ -> X |0⟩")
solve!(solver0)

println("\n|1⟩ -> X |1⟩")
solve!(solver1)

println("\nplotting results...")

# plot X gate on |0⟩ results

plot_wfn_ctrl_derivs(
    solver0,
    "plots/single_fluxonium_qubit/single_qstate_X_gate_on_0_basis_all.png",
    fig_title="X gate on 0 basis state",
)

plot_wfn(
    solver0,
    "plots/single_fluxonium_qubit/single_qstate_X_gate_on_0_basis_wfn.png",
    title="X gate on 0 basis state",
)

# plot X gate on |1⟩ results

plot_wfn_ctrl_derivs(
    solver1,
    "plots/single_fluxonium_qubit/single_qstate_X_gate_on_1_basis_all.png",
    fig_title="X gate on 1 basis state",
)

plot_wfn(
    solver1,
    "plots/single_fluxonium_qubit/single_qstate_X_gate_on_1_basis_wfn.png",
    title="X gate on 1 basis state",
)
