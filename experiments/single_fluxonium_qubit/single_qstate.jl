using QuantumControl

println("\npackages are loaded!")

# getting Hamiltonians from saved qutip objects
qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

N = 1001

# qubit basis states |0⟩ and |1⟩
ψ0 = [1, 0]
ψ1 = [0, 1]

gates = [:X, :Y, :Z, :H]

for gate in gates

    println("\nbeginning ALTRO solves for $gate gate...")

    prob_gate_on_0 = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψ0; # initial state
        N=N,
    )

    solver0 = ALTROSolver(prob_gate_on_0)

    prob_gate_on_1 = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψ1; # initial state
        N=N,
    )

    solver1 = ALTROSolver(prob_gate_on_1)


    println("\n|0⟩ -> $gate |0⟩")
    solve!(solver0)

    println("\n|1⟩ -> $gate |1⟩")
    solve!(solver1)

    println("\nplotting results...")

    # plot X gate on |0⟩ results

    plot_wfn_ctrl_derivs(
        solver0,
        "plots/single_fluxonium_qubit/single_qstate_$(gate)_gate_on_0_basis_all.png",
        fig_title="$gate gate on 0 basis state",
    )

    plot_wfn(
        solver0,
        "plots/single_fluxonium_qubit/single_qstate_$(gate)_gate_on_0_basis_wfn.png",
        title="$gate gate on 0 basis state",
    )

    # plot X gate on |1⟩ results

    plot_wfn_ctrl_derivs(
        solver1,
        "plots/single_fluxonium_qubit/single_qstate_$(gate)_gate_on_1_basis_all.png",
        fig_title="$gate gate on 1 basis state",
    )

    plot_wfn(
        solver1,
        "plots/single_fluxonium_qubit/single_qstate_$(gate)_gate_on_1_basis_wfn.png",
        title="$gate gate on 1 basis state",
    )
end
