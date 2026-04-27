using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using StaticArrays
using GLMakie
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
    return "ηᵣ = $(η_string) Pa s"
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
        fig = Figure(fontsize = 30, size = (800, 600) .* 2)
        ax = Axis(fig[1, 1], title = "Visco-elasto-viscoplastic model", xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]")
        if darkmode
            ax.xgridvisible = true
            ax.ygridvisible = true
            ax.xgridcolor = RGBAf(1, 1, 1, 1)
            ax.ygridcolor = RGBAf(1, 1, 1, 1)
            # ax.xgridwidth = 1
            # ax.ygridwidth = 2
            ax.spinewidth = 3
        end

        colors = darkmode ? (:tomato, :cyan, :gold, :lime, :magenta) : (:red, :blue, :orange, :green, :purple)
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

        yield_stress = first(results).yield_stress
        hlines!(
            ax,
            yield_stress / 1.0e6,
            color = darkmode ? :white : :black,
            linestyle = :dash,
            linewidth = 5,
            # label = "yield stress ($(round(yield_stress / 1.0e6, digits = 2)) MPa)",
            label = "yield stress)",
        )

        Legend(fig[1, 2], ax)
        colgap!(fig.layout, 20)
        display(fig)
    end

    dark_latex_theme = merge(
        theme_dark(),
        theme_latexfonts(),
        Theme(;
            textcolor = :white,
            Axis = (;
                titlecolor = :white,
                xlabelcolor = :white,
                ylabelcolor = :white,
                xticklabelcolor = :white,
                yticklabelcolor = :white,
                xgridvisible = true,
                ygridvisible = true,
                xtickcolor = :white,
                ytickcolor = :white,
                xgridcolor = RGBAf(1, 1, 1, 1),
                ygridcolor = RGBAf(1, 1, 1, 1),
                xgridwidth = 1,
                ygridwidth = 1,
                leftspinevisible = true,
                rightspinevisible = true,
                bottomspinevisible = true,
                topspinevisible = true,
                leftspinecolor = :white,
                rightspinecolor = :white,
                bottomspinecolor = :white,
                topspinecolor = :white,
                spinewidth = 3,
                xticksvisible = true,
                yticksvisible = true,
            ),
            Legend = (;
                labelcolor = :white,
                titlecolor = :white,
            ),
        ),
    )
    plot_theme = darkmode ? dark_latex_theme : theme_latexfonts()
    with_theme(() -> figure(darkmode = true), plot_theme)
end
