using GLMakie
using Printf

include("../rheologies/RheologyDefinitions.jl")

const SecYear = 365.25 * 24 * 3600

const εxx_pure_shear_3D = (1.0, -1.0, 0.0, 0.0, 0.0, 0.0)

second_invariant_3D(ε) = sqrt(0.5 * (ε[1]^2 + ε[2]^2 + ε[3]^2) + ε[4]^2 + ε[5]^2 + ε[6]^2)

function tensor_strain_rate_3D(εII; direction = εxx_pure_shear_3D)
    directionII = second_invariant_3D(direction)
    return @. εII * direction / directionII
end

vars_3D(εII, θ = 0.0; direction = εxx_pure_shear_3D) = (; ε = tensor_strain_rate_3D(εII; direction), θ)

function creep_rates(diffusion, dislocation, τ, env)
    ε̇_diff = compute_strain_rate(diffusion; τ = τ, T = env.T, P = env.P, f = env.f, d = env.d)
    ε̇_disl = compute_strain_rate(dislocation; τ = τ, T = env.T, P = env.P, f = env.f)
    return ε̇_diff, ε̇_disl
end

function creep_rate_derivative(r::DiffusionCreep, τ, env)
    iszero(τ) && return zero(τ)
    return r.n * compute_strain_rate(r; τ = τ, T = env.T, P = env.P, f = env.f, d = env.d) / τ
end

function creep_rate_derivative(r::DislocationCreep, τ, env)
    iszero(τ) && return zero(τ)
    return r.n * compute_strain_rate(r; τ = τ, T = env.T, P = env.P, f = env.f) / τ
end

function residual_and_jacobian(diffusion, dislocation, elastic, τ, εII_eff, env, dt)
    ε̇_diff, ε̇_disl = creep_rates(diffusion, dislocation, τ, env)
    residual = ε̇_diff + ε̇_disl + τ / (2 * elastic.G * dt) - εII_eff
    jacobian = creep_rate_derivative(diffusion, τ, env) +
               creep_rate_derivative(dislocation, τ, env) +
               1 / (2 * elastic.G * dt)
    return residual, jacobian
end

function solve_stress_newton(diffusion, dislocation, elastic, τ_guess, εII_eff, env, dt;
        atol = 1.0e-20, rtol = 1.0e-10, itermax = 50)

    τ = max(τ_guess, 0.0)
    r0, _ = residual_and_jacobian(diffusion, dislocation, elastic, τ, εII_eff, env, dt)
    tolerance = max(atol, rtol * abs(r0))

    for _ in 1:itermax
        r, drdτ = residual_and_jacobian(diffusion, dislocation, elastic, τ, εII_eff, env, dt)
        abs(r) ≤ tolerance && return τ

        Δτ = -r / drdτ
        α = 1.0
        while α > 1.0e-8
            τ_trial = max(τ + α * Δτ, 0.0)
            r_trial, _ = residual_and_jacobian(diffusion, dislocation, elastic, τ_trial, εII_eff, env, dt)
            abs(r_trial) ≤ (1 - 1.0e-4 * α) * abs(r) && break
            α *= 0.5
        end

        τ = max(τ + α * Δτ, 0.0)
    end

    error("Newton-Raphson stress solve did not converge after $itermax iterations; τ = $τ Pa")
end

function simulate_crustal_creep_newton_3D(diffusion, dislocation, elastic, vars, env;
        ntime = 501, dt = 250 * SecYear)

    time = zeros(ntime)
    τII = zeros(ntime)
    ε̇_diff = zeros(ntime)
    ε̇_disl = zeros(ntime)

    εII = second_invariant_3D(vars.ε)
    τ_guess = 2 * elastic.G * dt * εII

    for i in 2:ntime
        εII_eff = εII + τII[i - 1] / (2 * elastic.G * dt)
        τII[i] = solve_stress_newton(diffusion, dislocation, elastic, τ_guess, εII_eff, env, dt)
        ε̇_diff[i], ε̇_disl[i] = creep_rates(diffusion, dislocation, τII[i], env)

        τ_guess = τII[i]
        time[i] = time[i - 1] + dt
    end

    return (; time, τII, ε̇_diff, ε̇_disl)
end

function plot_crustal_creep_rates_3D(history, diffusion, dislocation, elastic, env)
    time_kyr = history.time / SecYear / 1.0e3
    dt = history.time[2] - history.time[1]
    τII = history.τII
    τII_MPa = τII ./ 1.0e6
    ε̇_diff = similar(τII)
    ε̇_disl = similar(τII)

    for i in eachindex(τII)
        ε̇_diff[i], ε̇_disl[i] = creep_rates(diffusion, dislocation, τII[i], env)
    end

    disl_fraction = ε̇_disl ./ (ε̇_diff .+ ε̇_disl)
    η_diff = τII ./ (2 .* ε̇_diff)
    η_disl = τII ./ (2 .* ε̇_disl)
    η_elastic = fill(elastic.G * dt, length(τII))
    η_eff = 1 ./ (1 ./ η_diff + 1 ./ η_disl + 1 ./ η_elastic)

    fig = Figure(fontsize = 24, size = (1300, 850), backgroundcolor = :black)
    axτ = Axis(fig[1, 1],
        xlabel = rich("t", " [kyr]", font = :italic),
        ylabel = rich("τ", subscript("II"), " [MPa]", font = :italic),
        backgroundcolor = :black,
    )
    axε = Axis(fig[1, 2],
        xlabel = rich("τ", subscript("II"), " [MPa]", font = :italic),
        ylabel = rich("ε̇", subscript("II"), " [s", superscript("-1"), "]", font = :italic),
        yscale = log10,
        backgroundcolor = :black,
    )
    axη = Axis(fig[2, 1],
        xlabel = rich("τ", subscript("II"), " [MPa]", font = :italic),
        ylabel = rich("η", " [Pa s]", font = :italic),
        yscale = log10,
        backgroundcolor = :black,
    )
    axf = Axis(fig[2, 2],
        xlabel = rich("τ", subscript("II"), " [MPa]", font = :italic),
        ylabel = "dislocation fraction",
        backgroundcolor = :black,
    )

    for ax in (axτ, axε, axη, axf)
        ax.titlecolor = :white
        ax.xlabelcolor = :white
        ax.ylabelcolor = :white
        ax.xticklabelcolor = :white
        ax.yticklabelcolor = :white
        ax.xgridcolor = (:white, 0.18)
        ax.ygridcolor = (:white, 0.18)
        ax.xminorgridcolor = (:white, 0.08)
        ax.yminorgridcolor = (:white, 0.08)
        ax.leftspinecolor = :white
        ax.rightspinecolor = :white
        ax.topspinecolor = :white
        ax.bottomspinecolor = :white
    end

    lines!(axτ, time_kyr, history.τII / 1.0e6, color = :white, linewidth = 3)

    history_positive = (history.τII .> 0) .& (history.ε̇_diff .> 0) .& (history.ε̇_disl .> 0)
    lines!(axε, τII_MPa[history_positive], history.ε̇_diff[history_positive], color = :dodgerblue3, linewidth = 3, label = "diffusion")
    lines!(axε, τII_MPa[history_positive], history.ε̇_disl[history_positive], color = :firebrick3, linewidth = 3, label = "dislocation")
    axislegend(axε, position = :rb, backgroundcolor = (:black, 0.0), framecolor = (:white, 0.25), labelcolor = :white)

    η_positive = (τII .> 0) .& isfinite.(η_eff) .& isfinite.(η_diff) .& isfinite.(η_disl)
    lines!(axη, τII_MPa[η_positive], η_eff[η_positive], color = :white, linewidth = 3, label = "effective")
    lines!(axη, τII_MPa[η_positive], η_diff[η_positive], color = :dodgerblue3, linewidth = 3, label = "diffusion")
    lines!(axη, τII_MPa[η_positive], η_disl[η_positive], color = :firebrick3, linewidth = 3, label = "dislocation")
    lines!(axη, τII_MPa[η_positive], η_elastic[η_positive], color = :gold, linewidth = 3, label = "elastic")
    axislegend(axη, position = :rt, backgroundcolor = (:black, 0.0), framecolor = (:white, 0.25), labelcolor = :white)

    fraction_positive = (τII .> 0) .& isfinite.(disl_fraction)
    lines!(axf, τII_MPa[fraction_positive], disl_fraction[fraction_positive], color = :limegreen, linewidth = 3)
    ylims!(axf, 0, 1)
    linkxaxes!(axη, axf)
    colgap!(fig.layout, 45)
    rowgap!(fig.layout, 35)

    display(fig)
    return fig
end

function setup_crustal_creep_newton_3D()
    R = 8.314462618

    env = (; T = 273.15 + 400, P = 500e6, f = 1.0, d = 1e-3)

    diffusion = DiffusionCreep(
        1,
        0.0,
        3.0,
        3e-21,
        160e3,
        8e-6,
        R,
    )

    dislocation = DislocationCreep(
        3.0,
        0.0,
        7e-26,
        190e3,
        8e-6,
        R,
    )

    elastic = IncompressibleElasticity(60e9)
    vars = vars_3D(1e-14, 0.0)

    return (; diffusion, dislocation, elastic, vars, env)
end

function main()
    setup = setup_crustal_creep_newton_3D()
    (; diffusion, dislocation, elastic, vars, env) = setup

    history = simulate_crustal_creep_newton_3D(diffusion, dislocation, elastic, vars, env)
    fig = plot_crustal_creep_rates_3D(history, diffusion, dislocation, elastic, env)

    println("Elastic + diffusion + dislocation creep, crust-like conditions, 3D stress tensor")
    println("Stress is solved with a scalar Newton-Raphson method.")
    println("T = $(env.T) K, P = $(env.P / 1e6) MPa, d = $(env.d * 1e3) mm, εII = $(second_invariant_3D(vars.ε)) s^-1")
    println("time [kyr]   τII [MPa]   diffusion [s^-1]   dislocation [s^-1]   creep fraction disl.")

    for i in 1:20:length(history.time)
        total_creep = history.ε̇_diff[i] + history.ε̇_disl[i]
        disl_fraction = iszero(total_creep) ? 0.0 : history.ε̇_disl[i] / total_creep
        @printf(
            "%10.1f%12.2f%20.4e%22.4e%21.3f\n",
            history.time[i] / SecYear / 1e3,
            history.τII[i] / 1e6,
            history.ε̇_diff[i],
            history.ε̇_disl[i],
            disl_fraction,
        )
    end

    return (; setup..., history, fig)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
