using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic, second_invariant

using GLMakie
import Statistics: mean

include("../rheologies/RheologyDefinitions.jl")

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))

function stress_time(c, ε, τ0, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    t_v  = zeros(ntime)
    
    P_e  = (0.0,)
    t    = 0.0

    for i in 2:ntime
        others   = (; dt = dt, P = others.P, P0 = P_e) # other non-differentiable variables needed to evaluate the state functions

        ε_corr = effective_strain_rate_correction(c, ε, τ0, others)
        ε_eff  = ε .+ ε_corr
        εII    = second_invariant(ε_eff...)

        vars   = merge(vars, (; ε = εII)) 

        x      = solve(c, x, vars, others, verbose = false)
        τII    = x[1] 
        τ0_vec = @. ε_eff * τII / εII
        τ0     = (τ0_vec[1],), (τ0_vec[2],),(τ0_vec[3],)

        τ1[i]  = τII
        
        t      += others.dt
        τ_an[i] = analytical_solution(second_invariant(ε...), t, c.leafs[2].G, c.leafs[1].η)
    
        t_v[i] = t
    end

    return t_v, τ1, τ_an
end

c, x, xnorm, ε, τ0, vars, args, others = let

    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)
    plastic = DruckerPrager(10e6, 30, 0)

    # Maxwell viscoelastic model
    # elastic --- viscous --- plastic
    c      = SeriesModel(viscous, elastic, plastic)

    ε      = 1e-14, -1e-14, 0e0
    τ0     = ((0e0, ), (0e0, ), (0e0,))
    vars   = (; ε = 1e-14, θ = 0e0)           # input variables (constant)
    args   = (; τ = 2e3, λ = 0e0)       # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e8, P = 1.0e6, P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x       = initial_guess_x(c, vars, args, others)
    char_τ  = plastic.C
    char_ε  = vars.ε 
    xnorm   = normalisation_x(c, char_τ, char_ε)

    c, x, xnorm, ε, τ0, vars, args, others
end

let
    t_v, τ, τ_an = stress_time(c,  ε, τ0, vars, x, xnorm, others; ntime = 1_500, dt = 1e8)

    function figure()
        SecYear = 3600 * 24 * 365.25
        fig = Figure(fontsize = 30, size = (800, 600) .* 2)
        ax  = Axis(fig[1, 1], title = "Visco-elasto-plastic model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")

        lines!(ax, t_v / SecYear / 1.0e3, τ_an / 1.0e6, color=:black, label = "viscoelastic analytical")
        scatter!(ax, t_v / SecYear / 1.0e3, τ / 1.0e6,  color=:red, label = "numerical")

        axislegend(ax, position = :rb)
        ax.xlabel = L"$t$ [kyr]"
        ax.ylabel = L"$\tau$ [MPa]"
        display(fig)
    end
    with_theme(figure, theme_latexfonts())
end
