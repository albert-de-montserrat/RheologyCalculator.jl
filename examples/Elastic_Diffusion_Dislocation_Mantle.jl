using RheologyCalculator
import RheologyCalculator: compute_strain_rate

using GLMakie
using Printf

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

const SecYear = 365.25 * 24 * 3600

function creep_rates(diffusion, dislocation, τ, env)
    ε̇_diff = compute_strain_rate(diffusion; τ = τ, T = env.T, P = env.P, f = env.f, d = env.d)
    ε̇_disl = compute_strain_rate(dislocation; τ = τ, T = env.T, P = env.P, f = env.f)
    return ε̇_diff, ε̇_disl
end

function simulate_mantle_creep(c, diffusion, dislocation, vars, x, xnorm, env;
        ntime = 501, dt = 250 * SecYear)

    time = zeros(ntime)
    τII = zeros(ntime)
    ε̇_diff = zeros(ntime)
    ε̇_disl = zeros(ntime)
    τ_elastic = (zero_stress_tensor_2D(),)

    for i in 2:ntime
        others = (; env..., dt = dt, τ0 = τ_elastic, P0 = (0.0,))
        x = solve(c, x, vars, others, xnorm0 = xnorm, verbose = false)

        τII[i] = x[1]
        ε̇_diff[i], ε̇_disl[i] = creep_rates(diffusion, dislocation, τII[i], env)
        τ_elastic = elastic_stress_history_2D(c, τII[i], vars.ε, τ_elastic, others)
        time[i] = time[i - 1] + dt
    end

    return (; time, τII, ε̇_diff, ε̇_disl, x)
end

function plot_mantle_creep_history(history, diffusion, dislocation, elastic, env)
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

c, diffusion, dislocation, elastic, vars, x, xnorm, env = let
    R = 8.314462618

    # Mantle-like conditions: upper-mantle temperature and pressure.
    env = (; T = 273.15 + 2500.0, P = 10.0e9, f = 1.0, d = 1.0e-3)

    # Dry olivine flow-law parameters after Hirth & Kohlstedt-style compilations.
    diffusion = DiffusionCreep(
        1,          # n, diffusion creep
        0.0,        # water fugacity exponent
        3.0,        # grain-size exponent
        1.58e-15,   # 10^9.2 MPa^-1 μm^3 s^-1 converted to Pa^-1 m^3 s^-1
        375e3,      # activation energy [J mol^-1]
        5e-6,       # activation volume [m^3 mol^-1]
        R,
    )

    dislocation = DislocationCreep(
        3.5,        # n, dislocation creep
        0.0,        # water fugacity exponent
        1.1e-16,    # 1.1e5 MPa^-3.5 s^-1 converted to Pa^-3.5 s^-1
        530e3,      # activation energy [J mol^-1]
        14e-6,      # activation volume [m^3 mol^-1]
        R,
    )

    elastic = IncompressibleElasticity(60e9)
    c = SeriesModel(diffusion, dislocation, elastic)

    vars = vars_2D(1e-15, 0.0)
    args = (; τ = 100e6, P = env.P)
    others = (; env..., dt = 250 * SecYear, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))

    x = initial_guess_x(c, vars, args, others)
    xnorm = normalisation_x(c, 150e6, second_invariant_2D(vars.ε))

    c, diffusion, dislocation, elastic, vars, x, xnorm, env
end

history = simulate_mantle_creep(c, diffusion, dislocation, vars, x, xnorm, env)
fig = plot_mantle_creep_history(history, diffusion, dislocation, elastic, env)
# save("mantle_history.png", fig)

println("Elastic + diffusion + dislocation creep, mantle-like conditions")
println("T = $(env.T) K, P = $(env.P / 1e9) GPa, d = $(env.d * 1e3) mm, εII = $(second_invariant_2D(vars.ε)) s^-1")
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
