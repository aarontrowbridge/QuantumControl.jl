using QuantumControl

println("\nloading saved Hamilitonian...")

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

# qubit basis states |0⟩ and |1⟩
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃1 = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

N = 1001
dt = 0.01

# X gate on two basis states |0⟩ and |1⟩
nqstates = 2
ψ̃i = [ψ̃0; ψ̃1]
ψ̃f = [ψ̃1; ψ̃0]

twoqstateprob = SingleQubitProblem(
    H_drift,
    H_drive,
    nqstates,
    ψ̃i,
    ψ̃f;
    dt=dt,
    N=N,
    ψ̃goal=false
)

println("\nsolving two qstate prob...")

solver = ALTROSolver(twoqstateprob; verbose=2)

solve!(solver)

println("\nplotting results...")

plot_two_wfns(
    solver,
    "plots/single_fluxonium_qubit/two_qstates_raw_X_gate_everything.png";
    show_dec_var=true,
    fig_title="X gate on basis states",
)
