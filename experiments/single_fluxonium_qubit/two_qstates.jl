using QuantumControl
using Altro
using TrajectoryOptimization
using StaticArrays

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
ψ̃1 = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

N = 1501
dt = 0.01

initprob0 = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃0, ψ̃1;
    N=N,
    dt=dt
    # U0_multiplier=0.1,
)

initprob1 = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃1, ψ̃0;
    N=N,
    dt=dt
    # U0_multiplier=0.1,
)

initsolver = ALTROSolver(initprob1; verbose=2)

solve!(initsolver)

plot_wfn_ctrl_derivs(initsolver, "plots/single_fluxonium_qubit/two_qstates_init_X_gate_on_1.png")

U0 = controls(initsolver)

println("\nsolving two quantum state problem (with U0 from init problem)...")

# X gate on two basis states |0⟩ and |1⟩
nqstates = 2
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f = @SVector [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0] # final quantum state isomorphism


twoqstateprob = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃0, ψ̃f;
    N=N,
    U0=U0,
)

dmodel, integrator = twoqstateprob.model

x0 = twoqstateprob.x0

t = range(0, N*dt, N)

xf = simulate(dmodel, x0, U0, t, dt, N)

plot_two_wfns(xf, t, "plots/single_fluxonium_qubit/two_qstates_X_gate_on_basis_qstates_simulated_from_1_solve.png"; show_dec=true, U=U0)



solver = ALTROSolver(twoqstateprob; verbose=2)

solve!(solver)

println("\nplotting results...")


plot_two_wfns(solver, "plots/single_fluxonium_qubit/two_qstates_wfns_ctrl.png"; show_dec=true)
