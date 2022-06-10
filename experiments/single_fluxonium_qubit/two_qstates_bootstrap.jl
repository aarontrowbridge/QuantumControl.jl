using QuantumControl

# plotting directory
plot_dir = "plots/single_fluxonium_qubit/"

println("\nloading saved Hamilitonian...")

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

println("\nsetting up single qstate prob...")

# qubit basis states |0⟩ and |1⟩
ψ0 = [1, 0]
ψ1 = [0, 1]

bootstrap_init_basis = ARGS[1]

gate = Symbol(ARGS[2])

bootstrap_nqstates = 1

T = 10.01
N = 1001
dt = T / N

if bootstrap_init_basis == "0"
    initprob = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψ0;
        N=N,
        dt=dt
    )
elseif bootstrap_init_basis == "1"
    initprob = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψ1;
        N=N,
        dt=dt
    )
else
    error("bootstrap_init_basis must be 0 or 1")
end

println("\nsolving init prob from $(bootstrap_init_basis) basis...")

initsolver = ALTROSolver(initprob; verbose=2)

solve!(initsolver)

plot_wfn_ctrl_derivs(
    initsolver,
    plot_dir*"two_qstates_bootstrap_$(gate)_gate_from_$(bootstrap_init_basis)_basis_init_all.png",
    fig_title="$gate gate on $bootstrap_init_basis basis state"
)

U0 = controls(initsolver)

# X gate on two basis states |0⟩ and |1⟩
nqstates = 2

ψ0s = [ψ0, ψ1, (ψ0 + im*ψ1)/√2, (ψ0 - ψ1)/√2]

twoqstateprob = SingleQubitProblem(
    H_drift,
    H_drive,
    gate,
    ψ0s;
    N=N,
    dt=dt,
    U0=U0,
    ψ̃goal=false
)

println("\nplotting 2 qstates with controls from init problem...")

X = simulate(twoqstateprob, U0)
t = range(0, N*dt, N)

plot_two_wfns(
    X,
    t,
    plot_dir*"two_qstates_bootstrap_$(gate)_gate_from_$(bootstrap_init_basis)_basis_sim_wfns.png";
    show_dec_var=true,
    U=U0,
    fig_title="$gate gate on basis states simulated from $bootstrap_init_basis basis solve"
)

println("\nsolving bootstrapped two qstate prob...")

solver = ALTROSolver(twoqstateprob; verbose=2)

solve!(solver)

println("\nplotting results...")

plot_two_wfns(
    solver,
    plot_dir*"two_qstates_bootstrap_$(gate)_gate_from_$(bootstrap_init_basis)_basis_result_all.png";
    show_dec_var=true,
    fig_title="bootstrapped $gate gate on basis states from $bootstrap_init_basis basis solve",
)
