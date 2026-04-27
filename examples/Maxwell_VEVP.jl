using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using StaticArrays
using GLMakie

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solutio vector
    τ1   = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (zero_stress_tensor_2D(),)
    P_e = (0.0,)
    t = 0.0
    for i in 2:ntime
        others   = (; dt = dt, τ0 = τ_e, P0 = P_e)      
        x        = solve(c, x, vars, others, verbose = true, xnorm0=xnorm)
        τ1[i]    = x[1]
        t       += others.dt
        τ_e      = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
    
        t_v[i] = t
    end
    return t_v, τ1
end

c, x, xnorm, vars, args, others, yield_stress = let
    viscous     = LinearViscosity(1e22)
    viscous_reg = LinearViscosity(1e20)
    # elastic     = Elasticity(10e9, 20e9)
    elastic     = IncompressibleElasticity(10e9)
    plastic     = DruckerPrager(10e6, 5, 0)

    # Maxwell visco-elasto-(visco-plastic) model
    p  = ParallelModel(plastic, viscous_reg)
    c  = SeriesModel(viscous, elastic, p)

    # input variables (constant)
    εᵢⱼ    = tensor_strain_rate_2D(1.0e-14)
    τ0ᵢⱼ   = (zero_stress_tensor_2D(),)
    vars   = (; ε = εᵢⱼ, θ = 0e0)
    # guess variables (we solve for these, differentiable)
    args   = (; τ = 2.0e3, λ = 0,  P = 1.0e6)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e8, τ0 = τ0ᵢⱼ, P0 = (0.0, ))

    x      = initial_guess_x(c, vars, args, others)
    char_τ = 1e6
    char_ε = second_invariant_2D(vars.ε)
    xnorm  = normalisation_x(c, char_τ, char_ε)

    yield_stress = args.P * plastic.sinϕ + plastic.C * plastic.cosϕ

    c, x, xnorm, vars, args, others, yield_stress
end

# eqs = RheologyCalculator.generate_equations(c)
# for eq in eqs
#     println(
#       "
#       self = $(eq.self)
#       fn   = $(eq.fn)
#       "
#     )
# end

let
    t_v, τ = stress_time(c, vars, x, xnorm, others; ntime = 1_500, dt = 1e8)
    darkmode = true

    function figure(; darkmode = false)
        SecYear = 3600 * 24 * 365.25
        fig = Figure(fontsize = 30, size = (800, 600) .* 2)
        ax  = Axis(fig[1, 1], title = "Visco-elasto-viscoplastic model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")
        if darkmode
            ax.xgridvisible = true
            ax.ygridvisible = true
            ax.xgridcolor = RGBAf(1, 1, 1, 1)
            ax.ygridcolor = RGBAf(1, 1, 1, 1)
            ax.xgridwidth = 8
            ax.ygridwidth = 8
        end

        linecolor = darkmode ? :tomato : :red
        lines!(ax, t_v / SecYear / 1.0e3, τ / 1.0e6, color = linecolor, label = "numerical", linewidth=5)
        hlines!(ax, yield_stress / 1.0e6, color = darkmode ? :cyan : :black, linestyle = :dash, linewidth = 5, label = "yield stress ($(round(yield_stress / 1.0e6, digits = 2)) MPa)")

        axislegend(ax, position = :rb)
        ax.xlabel = L"$t$ [kyr]"
        ax.ylabel = L"$\tau$ [MPa]"
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
    with_theme(() -> figure(darkmode = darkmode), plot_theme)
end
