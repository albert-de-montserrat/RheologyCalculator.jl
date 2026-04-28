using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using StaticArrays
using GLMakie
using LaTeXStrings

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
    elastic     = IncompressibleElasticity(30e9)
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
        fig = Figure(fontsize = 24, size = (1300, 850), backgroundcolor = darkmode ? :black : :white)
        ax  = Axis(fig[1, 1],
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

        linecolor = darkmode ? :dodgerblue3 : :red
        lines!(ax, t_v / SecYear / 1.0e3, τ / 1.0e6, color = linecolor, label = "numerical", linewidth = 3)
        hlines!(ax, yield_stress / 1.0e6, color = darkmode ? :white : :black, linestyle = :dash, linewidth = 3, label = "yield stress ($(round(yield_stress / 1.0e6, digits = 2)) MPa)")

        axislegend(ax, position = :rb, backgroundcolor = darkmode ? (:black, 0.0) : :white, framecolor = darkmode ? (:white, 0.25) : :black, labelcolor = darkmode ? :white : :black)
        display(fig)
    end
    figure(darkmode = darkmode)
end
