using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie
using LaTeXStrings

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8, verbose = false)
    τ1 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (zero_stress_tensor_2D(),)
    P_e = (0.0,)
    t = 0.0

    for i in 2:ntime
        others = (; others..., dt = dt, τ0 = τ_e, P0 = P_e)
        x = solve(c, x, vars, others, verbose = verbose, xnorm0 = xnorm, atol = 1.0e-10, rtol = 1e-10)
        τ1[i] = x[1]
        t += others.dt
        τ_e = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
        t_v[i] = t
    end

    return t_v, τ1
end

c1, x1, xnorm1, c2, x2, xnorm2, vars, others, yield_stress = let
    R = 8.314462618
    viscous = DislocationCreep(
        3.2,       # n, mafic granulite-style dislocation creep
        0.0,       # water fugacity exponent
        6.31e-22,  # Pa^-n s^-1
        244e3,     # activation energy [J mol^-1]
        4e-6,       # activation volume [m^3 mol^-1]
        R,
    )
    viscous_reg = LinearViscosity(1.0e20)
    elastic = IncompressibleElasticity(30.0e9)
    plastic = DruckerPrager(10.0e6, 5.0, 0.0)
    # plastic = DruckerPrager(Inf, 5.0, 0.0)

    # Maxwell visco-elasto-(visco-plastic) model with dislocation-creep viscous branch.
    p = ParallelModel(plastic, viscous_reg)
    c1 = SeriesModel(viscous, elastic, p)
    c2 = SeriesModel(viscous, elastic)

    εᵢⱼ = tensor_strain_rate_2D(1.0e-14)
    τ0ᵢⱼ = (zero_stress_tensor_2D(),)
    vars = (; ε = εᵢⱼ, θ = 0.0)
    args1 = (; τ = 2.0e3, λ = 0.0, P = 1.0e6)
    args2 = (; τ = 2.0e3)
    env = (; T = 273.15 + 400.0, P = 500.0e6, f = 1.0)
    env = (; T = 273.15 + 400.0, P = 500.0e6, f = 1.0)
    others = (; env..., dt = 1.0e8, τ0 = τ0ᵢⱼ, P0 = (0.0,))

    char_τ = 1.0e6
    char_ε = second_invariant_2D(vars.ε)
    x1 = initial_guess_x(c1, vars, args1, others)
    x2 = initial_guess_x(c2, vars, args2, others)
    xnorm1 = normalisation_x(c1, char_τ, char_ε)
    xnorm2 = normalisation_x(c2, char_τ, char_ε)
    yield_stress = env.P * plastic.sinϕ + plastic.C * plastic.cosϕ

    c1, x1, xnorm1, c2, x2, xnorm2, vars, others, yield_stress
end

let
    t_v1, τ1 = stress_time(c1, vars, x1, xnorm1, others; ntime = 8_000, dt = 1.0e8)
    t_v2, τ2 = stress_time(c2, vars, x2, xnorm2, others; ntime = 8_000, dt = 1.0e8)
    darkmode = true

    function figure(; darkmode = false)
        SecYear = 3600 * 24 * 365.25
        temperature_C = others.T - 273.15
        pressure_GPa = others.P / 1.0e9
        fig = Figure(fontsize = 24, size = (1300, 850), backgroundcolor = darkmode ? :black : :white)
        ax = Axis(fig[1, 1],
            title = "T = $(round(temperature_C; digits = 0)) °C, P = $(round(pressure_GPa; digits = 2)) GPa",
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

        lines!(ax, t_v1 / SecYear / 1.0e3, τ1 / 1.0e6,
            color = :deepskyblue2,
            label = L"V-E-VP",
            linewidth = 3,
        )
        lines!(ax, t_v2 / SecYear / 1.0e3, τ2 / 1.0e6,
            color = :orange,
            label = L"V-E",
            linewidth = 3,
        )
        hlines!(ax, yield_stress / 1.0e6,
            color = darkmode ? :white : :black,
            linestyle = :dash,
            linewidth = 3,
            label = L"\mathrm{yield\ stress}",
        )

        Legend(fig[1, 2], ax;
            backgroundcolor = darkmode ? (:black, 0.0) : :white,
            framecolor = darkmode ? (:white, 0.25) : :black,
            labelcolor = darkmode ? :white : :black,
        )
        colgap!(fig.layout, 35)
        display(fig)
        save("VEVP_disl.png", fig)
        return fig
    end

    figure(darkmode = darkmode)
end
