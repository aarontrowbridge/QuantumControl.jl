module PlottingUtils

export plot_wfn_ctrl_derivs
export plot_wfn
export plot_two_wfns

using CairoMakie
using Altro
using TrajectoryOptimization

# plot wavefunction components and a, ȧ, ä
function plot_wfn_ctrl_derivs(solver, filename, idx=1; fig_title=nothing)
    X = states(solver)
    U = controls(solver)
    t = gettimes(solver)

    N = length(X)

    ψ̃ = zeros(4, N)
    a = zeros(N)
    ȧ = zeros(N)
    ä = zeros(N)

    for k = 1:N
        ψ̃[:,k] = X[k][1 + (idx - 1)*4:idx*4]
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

    if !isnothing(fig_title)
        Label(fig[0,:], fig_title; textsize=30)
    end

    save(filename, fig)
end

# plot wavefunction components
function plot_wfn(solver, filename, idx=1; title=nothing)
    X = states(solver)
    t = gettimes(solver)

    N = length(X)

    ψ̃ = zeros(4, N)

    for k = 1:N
        ψ̃[:,k] = X[k][1 + (idx - 1)*4:idx*4]
    end

    fig = Figure(resolution=(600, 500))

    if isnothing(title)
        title="wavefunction components"
    else
        title=title
    end

    ax = Axis(fig[1, 1]; title=title, xlabel=L"t")

    series!(ax, t, ψ̃, labels=[L"\mathrm{Re} (\psi_0)",
                              L"\mathrm{Re} (\psi_1)",
                              L"\mathrm{Im} (\psi_0)",
                              L"\mathrm{Im} (\psi_1)"])
    axislegend(ax; position=:cb)

    save(filename, fig)
end

# plot two wavefunctions

function plot_two_wfns(solver::ALTROSolver, args...; show_dec_var=false, kwargs...)
    X = states(solver)
    t = gettimes(solver)
    U = nothing
    if show_dec_var
        U = controls(solver)
    end
    plot_two_wfns(X, t, args...; U=U, show_dec_var=show_dec_var, kwargs...)
end

function plot_two_wfns(X, t, filename, idxs=[1,2]; U=nothing, show_dec_var=false, fig_title=nothing)
    N = length(X)

    ψ̃₁ = zeros(4, N)
    ψ̃₂ = zeros(4, N)

    for k = 1:N
        ψ̃₁[:,k] = X[k][1 + (idxs[1] - 1)*4:idxs[1]*4]
        ψ̃₂[:,k] = X[k][1 + (idxs[2] - 1)*4:idxs[2]*4]
    end

    if !isnothing(U)
        fig = Figure(resolution=(1200, 1000))
    else
        fig = Figure(resolution=(1200, 500))
    end

    ax1 = Axis(fig[1, 1]; title="1st wavefunction components", xlabel=L"t")
    ax2 = Axis(fig[1, 2]; title="2nd wavefunction components", xlabel=L"t")

    series!(ax1, t, ψ̃₁; labels=[L"\mathrm{Re} (\psi_0)",
                                L"\mathrm{Re} (\psi_1)",
                                L"\mathrm{Im} (\psi_0)",
                                L"\mathrm{Im} (\psi_1)"])
    axislegend(ax1; position=:cb)

    series!(ax2, t, ψ̃₂; labels=[L"\mathrm{Re} (\psi_0)",
                                L"\mathrm{Re} (\psi_1)",
                                L"\mathrm{Im} (\psi_0)",
                                L"\mathrm{Im} (\psi_1)"])
    axislegend(ax2; position=:cb)

    if !isnothing(U) && show_dec_var
        ä = zeros(N)

        for k = 1:N-1
            ä[k+1] = U[k][1]
        end

        ax3 = Axis(fig[2, :]; title="decision variable", xlabel=L"t")

        lines!(ax3, t, ä; label=L"\frac{d^2a}{dt^2}(t)")
        axislegend(ax3; position=:rt)
    end

    if !isnothing(fig_title)
        Label(fig[0,:], fig_title; textsize=30)
    end

    save(filename, fig)
end



end
