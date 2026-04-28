using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using StaticArrays
using GLMakie
using LaTeXStrings
using Printf

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8, verbose = false)
    τ1 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (zero_stress_tensor_2D(),)
    P_e = (0.0,)
    t = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)
        x = solve(c, x, vars, others, verbose = verbose, xnorm0 = xnorm)
        τ1[i] = x[1]
        t += others.dt
        τ_e = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)

        t_v[i] = t
    end
    return t_v, τ1
end

function setup_model(η_reg)
    viscous = LinearViscosity(1e22)
    viscous_reg = LinearViscosity(η_reg)
    elastic = IncompressibleElasticity(10e9)
    plastic = DruckerPrager(10e6, 5, 0)

    p = ParallelModel(plastic, viscous_reg)
    c = SeriesModel(viscous, elastic, p)

    εᵢⱼ = tensor_strain_rate_2D(1.0e-14)
    τ0ᵢⱼ = (zero_stress_tensor_2D(),)
    vars = (; ε = εᵢⱼ, θ = 0e0)
    args = (; τ = 2.0e3, λ = 0, P = 1.0e6)
    others = (; dt = 1.0e8, τ0 = τ0ᵢⱼ, P0 = (0.0,))

    x = initial_guess_x(c, vars, args, others)
    char_τ = 1e6
    char_ε = second_invariant_2D(vars.ε)
    xnorm = normalisation_x(c, char_τ, char_ε)
    yield_stress = args.P * plastic.sinϕ + plastic.C * plastic.cosϕ

    return c, x, xnorm, vars, args, others, yield_stress
end

function eta_label(η)
    η_string = replace(@sprintf("%.0e", η), "e+" => "e")
    return latexstring("\\eta_r = $(η_string)\\ \\mathrm{Pa\\ s}")
end

let
    η_regs = (1e18, 1e19, 5e19, 1e20)
    darkmode = true
    ntime = 1_500
    dt = 1e8

    results = map(η_regs) do η_reg
        c, x, xnorm, vars, args, others, yield_stress = setup_model(η_reg)
        t_v, τ = stress_time(c, vars, x, xnorm, others; ntime = ntime, dt = dt)
        (; η_reg, t_v, τ, yield_stress)
    end

    function figure(; darkmode = false)
        SecYear = 3600 * 24 * 365.25
        fig = Figure(fontsize = 24, size = (1300, 850), backgroundcolor = darkmode ? :black : :white)
        ax = Axis(fig[1, 1],
            title = "Visco-elasto-viscoplastic model",
            xlabel = L"t\ [\mathrm{kyr}]",
            ylabel = L"\tau\ [\mathrm{MPa}]",
            backgroundcolor = darkmode ? :black : :white,
        )
        if darkmode
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

        yield_stress = first(results).yield_stress
        hlines!(
            ax,
            yield_stress / 1.0e6,
            color = darkmode ? (:white, 0.72) : :black,
            linestyle = :dash,
            linewidth = 3,
            label = L"\mathrm{yield\ stress}",
        )

        colors = darkmode ? ("#4CC9F0", "#F72585", "#F9C74F", "#90BE6D", "#B5179E") : (:red, :blue, :orange, :green, :purple)
        for (result, color) in zip(results, colors)
            lines!(
                ax,
                result.t_v / SecYear / 1.0e3,
                result.τ / 1.0e6,
                color = color,
                label = eta_label(result.η_reg),
                linewidth = 5,
            )
        end

        legend = Legend(fig[1, 2], ax, backgroundcolor = (:black, 0.0), framecolor = darkmode ? (:white, 0.25) : (:black, 0.25), labelcolor = darkmode ? :white : :black)
        colgap!(fig.layout, 35)
        display(fig)
        save("VEVP.png", fig)
    end

    figure(darkmode = darkmode)
end
