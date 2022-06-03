using QuantumControl

println("\nloading saved Hamilitonian...")

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

println("\nsetting up single qstate probs...")

# qubit basis states |0⟩ and |1⟩
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃1 = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

nqstates = 1

N = 1001
dt = 0.01

initprob0 = SingleQubitProblem(
    H_drift,
    H_drive,
    nqstates,
    ψ̃0,
    ψ̃1;
    N=N,
    dt=dt
    # U0_multiplier=0.1,
)

initprob1 = SingleQubitProblem(
    H_drift,
    H_drive,
    nqstates,
    ψ̃1,
    ψ̃0;
    N=N,
    dt=dt
    # U0_multiplier=0.1,
)

initprob = initprob0

println("\nsolving init prob...")

initsolver = ALTROSolver(initprob; verbose=2)

solve!(initsolver)

plot_wfn_ctrl_derivs(
    initsolver,
    "plots/single_fluxonium_qubit/two_qstates_init_X_gate_on_0.png",
    fig_title="X gate on 0 basis state"
)

U0 = controls(initsolver)

println("\nsolving two quantum state problem with a(t) from init problem...")

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
    N=N,
    U0=U0,
    ψ̃goal=false
)

println("\nplotting with controls from init problem...")

X = simulate(twoqstateprob, U0)
t = range(0, N*dt, N)

plot_two_wfns(
    X,
    t,
    "plots/single_fluxonium_qubit/two_qstates_X_gate_on_basis_qstates_simulated_from_0_solve.png";
    show_dec_var=true,
    U=U0,
    fig_title="X gate on basis states simulated from 1st solve"
)

println("\nsolving two qstate prob...")

solver = ALTROSolver(twoqstateprob; verbose=2)

solve!(solver)

println("\nplotting results...")

plot_two_wfns(
    solver,
    "plots/single_fluxonium_qubit/two_qstates_wfns_ctrl.png";
    show_dec_var=true,
    fig_title="X gate on basis states",
)
