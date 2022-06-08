using QuantumControl

import Base.Threads: @threads

println("\nloading saved Hamilitonian...")

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

T = 10.01
N = 1001
dt = T / N

gate = :X

# X gate on two basis states |0⟩ and |1⟩
# nqstates = 2
# ψ̃i = [ψ̃0; ψ̃1]
# ψ̃f = [ψ̃1; ψ̃0]

ψi = [[1, 0], [0, 1]]

gates = [:X, :Y, :Z, :H]

solvers = Vector{Altro.AbstractSolver}(undef, length(gates))

@threads for i = 1:length(gates)
    gate = gates[i]

    twoqstateprob = SingleQubitProblem(
        H_drift,
        H_drive,
        gate,
        ψi;
        dt=dt,
        N=N,
        # Qᵢᵢ=0.1,
        # ctrl_cost=0.1,
        # ctrl_cost_multiplier=100.0,
        # state_cost_multiplier=1000.0,
        # dec_cost_multiplier=10.0,

    )

    println("\nsolving two qstate $gate gate prob...")

    solver = ALTROSolver(
        twoqstateprob;
        max_cost_value=1e10,
        iterations=2000,
        iterations_outer=100
    )

    solve!(solver)
    solvers[i] = solver
end

println("\nplotting results...")

for (gate, solver) in zip(gates, solvers)
    plot_two_wfns(
        solver,
        "plots/single_fluxonium_qubit/two_qstates_raw_$(gate)_gate_wfns_dec.png";
        show_dec_var=true,
        fig_title="$gate gate on basis states",
    )
end
