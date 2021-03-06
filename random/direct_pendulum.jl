


# PREAMBLE

# PKG_SETUP

# ## Setup

using LinearAlgebra
using Plots
using DirectTrajectoryOptimization

# ##
eval_hess_lag = true

# ## horizon
T = 11

# ## acrobot
num_state = 2
num_action = 1
num_parameter = 0

function pendulum(x, u, w)
    mass = 1.0
    length_com = 0.5
    gravity = 9.81
    damping = 0.1

    [
        x[2],
        (u[1] / ((mass * length_com * length_com))
            - gravity * sin(x[1]) / length_com
            - damping * x[2] / (mass * length_com * length_com))
    ]
end

function midpoint_implicit(y, x, u, w)
    h = 0.05 # timestep
    y - (x + h * pendulum(0.5 * (x + y), u, w))
end

# ## model
dt = Dynamics(
    midpoint_implicit,
    num_state,
    num_state,
    num_action,
    eval_hess=eval_hess_lag)
dynamics = [dt for t = 1:T-1]

# ## initialization
x1 = [0.0; 0.0]
xT = [π; 0.0]

# ## objective
ot = (x, u, w) -> 0.1 * dot(x[1:2], x[1:2]) + 0.1 * dot(u, u)
oT = (x, u, w) -> 0.1 * dot(x[1:2], x[1:2])
ct = Cost(ot, num_state, num_action,
    eval_hess=eval_hess_lag)
cT = Cost(oT, num_state, 0,
    eval_hess=eval_hess_lag)
objective = [[ct for t = 1:T-1]..., cT]

# ## constraints
bnd1 = Bound(num_state, num_action)
bndt = Bound(num_state, num_action)
bndT = Bound(num_state, 0)
bounds = [bnd1, [bndt for t = 2:T-1]..., bndT]

# ## initial
con1 = Constraint((x, u, w) -> x - x1, num_state, num_action,
    eval_hess=eval_hess_lag)
conT = Constraint((x, u, w) -> x - xT, num_state, num_action,
    eval_hess=eval_hess_lag)
constraints = [con1, [Constraint() for t = 2:T-1]..., conT]

# ## problem
solver = Solver(dynamics, objective, constraints, bounds,
    eval_hess=eval_hess_lag)

# ## initialize
x_interpolation = linear_interpolation(x1, xT, T)
u_guess = [1.0 * randn(num_action) for t = 1:T-1]

initialize_states!(solver, x_interpolation)
initialize_controls!(solver, u_guess)

# ## solve
@time solve!(solver)

# ## solution
x_sol, u_sol = get_trajectory(solver)

@show x_sol[1]
@show x_sol[T]

# ## state
plot(hcat(x_sol...)')

# ## control
plot(hcat(u_sol[1:end-1]..., u_sol[end-1])',
    linetype=:steppost)
