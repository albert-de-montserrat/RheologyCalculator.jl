using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using Statistics: mean

using GLMakie
GLMakie.activate!(; visible = false)

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

# Analytical solution for the Mixed Kelvin-Voigt and Maxwell body
#
# Model: SeriesModel(η₁, ParallelModel(η₂, SeriesModel(η₃, G)))
#
#         ┌──── η₂ ─────┐
#         │             │
#  η₁ ─── │             │
#         │             │
#         └─ η₃ ─── G ──┘
#
# Under constant strain rate εII starting from zero stress, the global
# stress obeys a first-order linear ODE with the exact solution
#
#   τ(t) = τ∞ - (τ∞ - τ₀) exp(-t / t_relax)
#
# τ₀      = 2η₁η₂/(η₁+η₂) εII          (t=0: spring unloaded)
# τ∞      = 2η₁(η₂+η₃)/(η₁+η₂+η₃) εII  (long-time equilibrium)
# t_relax = (η₁+η₂)η₃ / (G(η₁+η₂+η₃))  (relaxation timescale)
function analytical_solution(t, εII, η1, η2, η3, G)
    τ_init  = 2η1 * η2 / (η1 + η2) * εII
    τ_inf   = 2η1 * (η2 + η3) / (η1 + η2 + η3) * εII
    t_relax = (η1 + η2) * η3 / (G * (η1 + η2 + η3))
    return τ_inf - (τ_inf - τ_init) * exp(-t / t_relax)
end

function stress_time(c, vars, x, xnorm; ntime = 200, dt = 1.0e9)
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    t_v  = zeros(ntime)
    τ_e  = (zero_stress_tensor_2D(),)
    P_e  = (0.0,)
    t    = 0.0
    εII  = second_invariant_2D(vars.ε)
    # access model parameters for the analytical solution
    η1 = c.leafs[1].η
    η2 = c.branches[1].leafs[1].η
    η3 = c.branches[1].branches[1].leafs[1].η
    G  = c.branches[1].branches[1].leafs[2].G
    for i in 2:ntime
        others  = (; dt = dt, τ0 = τ_e, P0 = P_e)
        x       = solve(c, x, vars, others; xnorm0 = xnorm)
        # The elastic element sits inside the inner SeriesModel; pass the full x
        # SVector so that compute_stress_elastic extracts the correct spring stress.
        τ_e     = elastic_stress_history_2D(c, x, vars.ε, τ_e, others)
        t      += dt
        τ1[i]   = x[1]
        τ_an[i] = analytical_solution(t, εII, η1, η2, η3, G)
        t_v[i]  = t
    end
    return t_v, τ1, τ_an
end

c, x, xnorm, vars, args, others = let
    η1 = LinearViscosity(1e22)
    η2 = LinearViscosity(1e21)
    η3 = LinearViscosity(1e21)
    el = IncompressibleElasticity(1e10)

    c = SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))

    εII    = 1.0e-14
    vars   = vars_2D(εII)
    args   = (; τ = 2.0e7, P = 0.0)
    others = (; dt = 1.0e6, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))
    # others = (; dt = 1.0e9, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))

    x = initial_guess_x(c, vars, args, others)

    # long-time equilibrium stress used as characteristic scale for normalisation
    τ_char = 2η1.η * (η2.η + η3.η) / (η1.η + η2.η + η3.η) * εII
    xnorm  = normalisation_x(c, τ_char, εII)

    c, x, xnorm, vars, args, others
end

let
    function figure()
        dt = 1e10 .* [1.0, 1/2, 1/4, 1/8, 1/64]
        nt = 100 .* [1.0, 2, 4, 8, 64]

        ϵ  = zero(dt)
        t_v = τ = τ_an = nothing

        for it in eachindex(dt)
            t_v, τ, τ_an = stress_time(c, vars, x, xnorm; ntime = Int64(nt[it]), dt = dt[it])
            ϵ[it] = (100 .* abs.(τ .- τ_an)./τ_an)[2:end] |> mean # average relative error in percent, ignoring the first point
        end

        θ      = log(ϵ[1]/ϵ[2]) / log(dt[1]/dt[2])
        dt_arr = LinRange(dt[1], dt[end], 100)
        ϵ_arr  = ϵ[1] .* (dt_arr ./ dt[1]).^θ

        SecYear = 3600 * 24 * 365.25
        fig = Figure(fontsize = 30, size = (800, 600))
        ax1 = Axis(fig[1, 1], title = L"$$Mixed KV–Maxwell model", xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]")
        lines!(ax1, t_v / SecYear / 1.0e3, τ_an / 1.0e6, color = :black, label = "analytical")
        scatter!(ax1, t_v[1:5:end] / SecYear / 1.0e3, τ[1:5:end] / 1.0e6, color = :red, label = "numerical")
        axislegend(ax1, position = :rb, labelsize = 18)

        ax2 = Axis(fig[2, 1], title = L"$$Convergence", xlabel = L"$\log_{10}$ $\frac{1}{dt}$ [1/s]", ylabel = L"$\log_{10}$ mean relative error [%]")
        lines!(ax2, log10.(1 ./ dt_arr), log10.(ϵ_arr), color = :black, label = "1st order")
        scatter!(ax2, log10.(1 ./ dt), log10.(ϵ), color = :black, label = "mean rel. error")
        axislegend(ax2, position = :rt, labelsize = 18)

        GLMakie.save("docs/assets/Maxwell_KV_Maxwell.png", fig)
        display(fig)
    end
    with_theme(figure, theme_latexfonts())
end
