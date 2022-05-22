# QuantumControl.jl

this package aims to provide an interface between the python package QuTiP and the Altro.jl trajectory optimization package.

the goal is to do multi-qubit quantum optimal control quickly and robustly. 

## installation

to use this package, clone it, and then, from a julia REPL in the cloned directory run

`(@v1.7) pkg> activate .`

`(QuantumControl) pkg> instantiate`

## usage

right now, all that is implemented is 

* functions to load saved qutip objects into julia
* a `RobotDynamics.ContinuousDynamics` subtype model for a multi-qubit system
* dynamics functions for this model which convert complex objects to isomporphic real objects

## quantum optimal control

in the future this section will discuss the theory of QOC. for now I will just demonstrate the new ability to render $\LaTeX$ in github markdown.

$$
i \hbar \partial_t \ket{\psi} = \hat H \ket{\psi}
$$
