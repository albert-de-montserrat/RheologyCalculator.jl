using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic, second_invariant

using GLMakie
import Statistics: mean

include("../rheologies/RheologyDefinitions.jl")

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))

function stress_time(c, ε, τ0, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    t_v  = zeros(ntime)
    # τ_e  = (0.0,)
    P_e  = (0.0,)
    t    = 0.0
    others = (; dt = dt, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

    x      = initial_guess_x(c, vars, args, others)

    for i in 2:ntime
        others = (; dt = dt, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

        ε_corr = effective_strain_rate_correction(c, ε, τ0, others)
        # @show  τ0
        ε_eff = ε .+ ε_corr
        εII   = second_invariant(ε_eff...)

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


c, x, ε, τ0, vars, args, others = let

    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c      = SeriesModel(viscous, elastic)

    ε      = 1e-14, -1e-14, 0e0
    τ0     = ((0e0, ), (0e0, ), (0e0,))
    vars   = (; ε = 1e-14, θ = 0e0)           # input variables (constant)
    args   = (; τ = 2e3, P = 1.0e0)       # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e9, P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, ε, τ0, vars, args, others
end

let
    function figure()
        dt = 1e10 .* [1.0, 1/2, 1/4, 1/8]
        nt = 1_000 .* [1.0, 2, 4, 8]
        ϵ  = zero(dt) 

        for it in eachindex(dt)
            t_v, τ, τ_an = stress_time(c, ε, τ0, vars, x; ntime = Int64(nt[it]), dt = dt[it])
            ϵ[it] = maximum(abs.(τ .- τ_an))

            # Order
            θ      = log(ϵ[1]/ϵ[2]) / log(dt[1]/dt[2])
            dt_arr = LinRange(dt[end], dt[1], 100)
            ϵ_arr  = 10 .^(log10.(ϵ[1]) .- θ.*( log10.(1 ./ (dt_arr)) .- log10.(1 ./ (dt_arr[end]))))

            SecYear = 3600 * 24 * 365.25
            fig = Figure(fontsize = 30, size = (800, 600) .* 1)

            ax1 = Axis(fig[1, 1], title = L"$$Maxwell model", xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]")
            lines!(ax1, t_v / SecYear / 1.0e3, τ_an / 1.0e6, color=:black, label = "analytical")
            scatter!(ax1, t_v[1:1000:end] / SecYear / 1.0e3, τ[1:1000:end] / 1.0e6,  color=:red, label = "numerical")
            axislegend(ax1, position = :rb, labelsize=18)
            
            ax2 = Axis(fig[2, 1], title = L"$$Convergence", xlabel = L"$\log_{10}$ $\frac{1}{dt} $ [1/s]", ylabel = L"$\log_{10}$ $ϵ$ [Pa]")
            lines!(ax2, log10.(1 ./ dt_arr), log10.(ϵ_arr), color=:black, label="1st order")
            scatter!(ax2, log10.(1 ./ dt), log10.(ϵ), color=:black, label="numerics")
            axislegend(labelsize=18)

            save("docs/assets/Maxwell_VE_model.png", fig)
            display(fig)
        end
    end
    with_theme(figure, theme_latexfonts())
end
