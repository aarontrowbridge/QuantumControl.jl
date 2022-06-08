using QuantumControl

import Base.Threads: @threads

# plotting directory
plot_dir = "plots/single_fluxonium_qubit/"

println("\npackages are loaded!")

# getting Hamiltonians from saved qutip objects
qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

# qubit basis states |0⟩ and |1⟩
ψ0 = [1, 0]
ψ1 = [0, 1]

# number of time steps
T = 10.01
N = 1001
dt = T / N

plot_wfn = false

äbound = true

fidelity_cost = true

# loop over all single qubit gates
gates = [:X, :Y, :Z, :H]
# gates = [:X]

solver0s = Vector{Altro.AbstractSolver}(undef, length(gates))
solver1s = Vector{Altro.AbstractSolver}(undef, length(gates))

@threads for i = 1:length(gates)
    gate = gates[i]

    println("\nbeginning ALTRO solves for $gate gate...")

    prob_gate_on_0 = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψ0; # initial state
        dt=dt,
        N=N,
        fidelity_cost=fidelity_cost,
        äbound=äbound,
    )

    solver0 = ALTROSolver(prob_gate_on_0)

    prob_gate_on_1 = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψ1; # initial state
        dt=dt,
        N=N,
        fidelity_cost=fidelity_cost,
        äbound=äbound,
    )

    solver1 = ALTROSolver(prob_gate_on_1)


    println("\n|0⟩ -> $gate |0⟩")
    solve!(solver0)
    solver0s[i] = solver0

    println("\n|1⟩ -> $gate |1⟩")
    solve!(solver1)
    solver1s[i] = solver1
end

println("\nplotting results...")

for (gate, solver0, solver1) in zip(gates, solver0s, solver1s)

    # plot gate on |0⟩ results

    plot_wfn_ctrl_derivs(
        solver0,
        plot_dir*"single_qstate_$(gate)_gate_on_0_basis_all_T_$(T)_$(fidelity_cost ? "fidelity" : "LQR").png",
        fig_title="$gate gate on 0 basis state",
    )

    # plot gate on |1⟩ results

    plot_wfn_ctrl_derivs(
        solver1,
        plot_dir*"single_qstate_$(gate)_gate_on_1_basis_all_T_$(T)_$(fidelity_cost ? "fidelity" : "LQR").png",
        fig_title="$gate gate on 1 basis state",
    )

    if plot_wfn
        plot_wfn(
            solver0,
            plot_dir*"single_qstate_$(gate)_gate_on_0_basis_wfn_T_$(T)_$(fidelity_cost ? "fidelity" : "LQR").png",
            title="$gate gate on 0 basis state",
        )

        plot_wfn(
            solver1,
            plot_dir*"single_qstate_$(gate)_gate_on_1_basis_wfn_T_$(T)_$(fidelity_cost ? "fidelity" : "LQR").png",
            title="$gate gate on 1 basis state",
        )
    end
end
