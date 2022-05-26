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

model = MultiQubitSystem(H_drift, H_drive, nqstates)

dmodel = RD.DiscretizedDynamics{RD.RK4}(model)

dt = 0.01 # time step
N = 101 # number of knot points
tf = dt * (N - 1) # final time

n, m = RD.dims(model)

nd = n - model.nstates

# Set initial and final states
q0 = @SVector [1.0,  0, 0,    0] # initial quaternion
qf = @SVector [1/√2, 0, 1/√2, 0] # final quaternion

# ∫a, a, da should start and stop at 0
x0 = [q0; @SVector zeros(nd)]
xf = [qf; @SVector zeros(nd)]

# random initial guess for control inputs
U0 = [@SVector randn(m) for k = 1:N-1]

#Set up quadratic objective function
Q = Diagonal(@SVector fill(1.0, n))
R = Diagonal(@SVector fill(0.1, m))
Qf = 100 * Q
obj = LQRObjective(Q, R, Qf, xf, N, checks=false);

# penalize control at first and last time step more
costfun = LQRCost(Q, 100 * R, xf)
obj.cost[1] = costfun
obj.cost[N-1] = costfun


# set up contraints
cons = ConstraintList(n, m, N)

# final goal
# goalcon = GoalConstraint(xf)
# add_constraint!(cons, goalcon, N)

# maintain normalized statevector
normcon = NormConstraint(n, m, 1.0, Equality(), :state)
add_constraint!(cons, normcon, N)

prob = Problem(dmodel, obj, x0, tf; constraints=cons, U0=U0, xf=xf)

solver = ALTROSolver(prob; verbose=true, projected_newton=true)

solve!(solver)

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
