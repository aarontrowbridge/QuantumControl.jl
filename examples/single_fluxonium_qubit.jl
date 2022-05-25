using QuantumControl
using Altro
using TrajectoryOptimization
using RobotDynamics
using LinearAlgebra
using StaticArrays
using Random

import TrajectoryOptimization as TO
import RobotDynamics as RD

qutip_obj_dir = "qutip_saved_objects/single_fluxonium_qubit/"

H_drift_path = qutip_obj_dir * "H_drift.qu"
H_drive_path = qutip_obj_dir * "H_drive.qu"

H_drift = get_qutip_matrix(H_drift_path)
H_drive = get_qutip_matrix(H_drive_path)

nqstates = 1

model = MultiQubitSystem(H_drift, H_drive, nqstates; isodynamics=false)

dmodel = RD.DiscretizedDynamics{RD.RK4}(model)

dt = 0.01 # time step
N = 101 # number of knot points
tf = dt * (N - 1) # final time

x0 = SA[1.0, 0, 0, 0, 0, 0, 0]
xf = SA[1/sqrt(2), 0, 1/sqrt(2), 0, 0, 0, 0]

n, m = RD.dims(model)

#Set up quadratic objective function
Q = 0.001 * Diagonal(@SVector ones(n))
R = 10.0 * Diagonal(@SVector ones(m))
Qf = 1.0 * Diagonal(@SVector ones(n))
obj = LQRObjective(Q,R,Qf,xf,N);

#Set up input constriants
goal = GoalConstraint(xf)
# con_goal = ConstraintVals(goal, N:N)
con_list = ConstraintList(n, m, N)
add_constraint!(con_list, goal, N)

prob = Problem(dmodel, obj, x0, tf; constraints=con_list)

solver = ALTROSolver(prob; verbose=true, projected_newton=true)

solve!(solver)

println("solved!")

X̃ = states(solver)
X = zeros(4,N)
U = zeros(N)
U̇ = zeros(N)
for k = 1:N
      X[:,k] = X̃[k][1:4]
      U[k] = X̃[k][6]
      U̇[k] = X̃[k][5]
end

using Plots

plot(X',xlabel="time step",title="State Trajectory",label=["Re(x1)" "Im(x1)" "Re(x2)" "Im(x2)"])
savefig("plots/single_fluxonium_qubit_state_trajectory.png")
plot(U,xlabel="time step",title="Control Trajectory",label="u")
savefig("plots/single_fluxonium_qubit_control_trajectory.png")
plot(U̇,xlabel="time step",title="Control Trajectory",label="udot")
savefig("plots/single_fluxonium_qubit_control_derivative_trajectory.png")
