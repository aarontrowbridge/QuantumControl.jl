using QuantumControl

# plotting directory
plot_dir = "plots/single_fluxonium_qubit/"

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

bootstrap_init_basis = parse(Int, ARGS[1])

bootstrap_nqstates = 1

N = 1001
dt = 0.01

initprob0 = SingleQubitProblem(
    H_drift,
    H_drive,
    bootstrap_nqstates,
    ψ̃0,
    ψ̃1;
    N=N,
    dt=dt
    # U0_multiplier=0.1,
)

initprob1 = SingleQubitProblem(
    H_drift,
    H_drive,
    bootstrap_nqstates,
    ψ̃1,
    ψ̃0;
    N=N,
    dt=dt
    # U0_multiplier=0.1,
)

if bootstrap_init_basis == 0
    initprob = initprob0
elseif bootstrap_init_basis == 1
    initprob = initprob1
else
    error("bootstrap_init_basis must be 0 or 1")
end

println("\nsolving init prob from $(bootstrap_init_basis) basis...")

initsolver = ALTROSolver(initprob; verbose=2)

solve!(initsolver)

plot_wfn_ctrl_derivs(
    initsolver,
    plot_dir*"two_qstates_bootstrap_from_$(bootstrap_init_basis)_basis__X_gate_init_all.png",
    fig_title="X gate on 0 basis state"
)

U0 = controls(initsolver)

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

println("\nplotting 2 qstates with controls from init problem...")

X = simulate(twoqstateprob, U0)
t = range(0, N*dt, N)

plot_two_wfns(
    X,
    t,
    plot_dir*"two_qstates_bootstrap_from_$(bootstrap_init_basis)_basis__X_gate_sim_wfns.png";
    show_dec_var=true,
    U=U0,
    fig_title="X gate on basis states simulated from 1st solve"
)

println("\nsolving bootstrapped two qstate prob...")

solver = ALTROSolver(twoqstateprob; verbose=2)

solve!(solver)

println("\nplotting results...")

plot_two_wfns(
    solver,
    plot_dir*"two_qstates_bootstrap_from_$(bootstrap_init_basis)_basis__X_gate_result_all.png";
    show_dec_var=true,
    fig_title="X gate on basis states",
)
