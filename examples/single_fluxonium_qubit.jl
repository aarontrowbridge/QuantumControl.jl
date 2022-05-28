using QuantumControl
using Altro
using TrajectoryOptimization
using RobotDynamics
using LinearAlgebra
using StaticArrays
using Random
using CairoMakie

import TrajectoryOptimization as TO
import RobotDynamics as RD

println("\npackages are loaded!")

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

nqstates = 1

# X gate
ψ̃0 = @SVector [1.0, 0.0, 0.0, 0.0] # initial quantum state isomorphism
ψ̃f = @SVector [0.0, 1.0, 0.0, 0.0] # final quantum state isomorphism

# # Y gate
# Yψ̃0 = @SVector [1.0, 0.0, 1.0, 1.0] # initial quantum state isomorphism

prob = SingleQubitProblem(H_drift, H_drive, nqstates, ψ̃0, ψ̃f; N=1001)

solver = ALTROSolver(prob; verbose=true, projected_newton=true)

println("\nbeginning ALTRO solve...")

solve!(solver)

println("\nplotting results...")

X = states(solver)
U = controls(solver)
t = gettimes(solver)

N = length(X)

ψ̃ = zeros(4, N)
a = zeros(N)
ȧ = zeros(N)
ä = zeros(N)

for k = 1:N
    ψ̃[:,k] = X[k][1:4]
    a[k] = X[k][6]
    ȧ[k] = X[k][7]
end

for k = 1:N-1
    ä[k+1] = U[k][1]
end

fig = Figure(resolution=(1200, 1000))

ax1 = Axis(fig[1, 1]; title="wavefunction components", xlabel=L"t")
ax2 = Axis(fig[1, 2]; title="control variable", xlabel=L"t")
ax3 = Axis(fig[2, 1]; title="first derivative of control", xlabel=L"t")
ax4 = Axis(fig[2, 2]; title="second derivative of control (decision variable)", xlabel=L"t")

series!(ax1, t, ψ̃; labels=[L"\mathrm{Re} (\psi_0)",
                           L"\mathrm{Re} (\psi_1)",
                           L"\mathrm{Im} (\psi_0)",
                           L"\mathrm{Im} (\psi_1)"])
axislegend(ax1; position=:cb)
lines!(ax2, t, a; label=L"a(t)")
axislegend(ax2; position=:rt)
lines!(ax3, t, ȧ; label=L"\frac{da}{dt}(t)")
axislegend(ax3; position=:rt)
lines!(ax4, t, ä; label=L"\frac{d^2a}{dt^2}(t)")
axislegend(ax4; position=:rt)

save("plots/single_fluxonium_qubit.png", fig)

fig_wfn = Figure(resolution=(600, 500))
ax_wfn = Axis(fig_wfn[1, 1]; title="X gate wavefunction components", xlabel=L"t")

series!(ax_wfn, t, ψ̃, labels=[L"\mathrm{Re} (\psi_0)",
                              L"\mathrm{Re} (\psi_1)",
                              L"\mathrm{Im} (\psi_0)",
                              L"\mathrm{Im} (\psi_1)"])
axislegend(ax_wfn; position=:cb)

save("plots/single_fluxonium_qubit_wfn.png", fig_wfn)
