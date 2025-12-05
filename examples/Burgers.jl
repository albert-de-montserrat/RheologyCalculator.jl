using ForwardDiff, RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("RheologyDefinitions.jl")

using GLMakie

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1 = zeros(ntime)
    τ2 = zeros(ntime)
    P1 = zeros(ntime)
    P2 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (0.0, 0.0)
    P_e = (0.0, 0.0)
    t = 0.0
    for i in 2:ntime
        #global τ_e, P_e, x, vars, others, t
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

        x = solve(c, x, vars, others)
        τ_e = compute_stress_elastic(c, x, others)
        P_e = compute_pressure_elastic(c, x, others)

        @inbounds τ1[i] = τ_e[1]
        @inbounds P1[i] = P_e[1]
        # if length(τ_e) > 1
            @inbounds τ2[i] = τ_e[2]
            @inbounds P2[i] = P_e[2]
        # end
        t += others.dt
        t_v[i] = t
    end
    return t_v, τ1, τ2, P1, P2, x
end

# Analytical solution for Burgers model
function simulate_series_Burgers_model(E1, η1, E2, η2, ε̇, N, dt)
    t    = LinRange(0., N*dt, N)
    σ    = zeros(N)          # Stress
    ε_KV = zeros(N)          # Strain in Kelvin–Voigt element
    σ_KV = zeros(N)
    σ_spring = zeros(N) # Stress in the spring of the Kelvin-Voigt element
    for i in 2:N
        # Previous values
        ε_KV_prev = ε_KV[i - 1]
        σ_prev = σ[i - 1]
        dεKVdt = (σ_prev - E2 * ε_KV_prev) / η2    # Kelvin–Voigt strain rate
        ε_KV[i] = ε_KV_prev + dt * dεKVdt           # Update ε_KV
        dσdt = E1 * (ε̇ - σ_prev / η1 - dεKVdt)   # Stress rate from Maxwell element
        σ[i] = σ_prev + dt * dσdt                # Update stress
        σ_KV[i] = E2 * ε_KV[i] + η2 * dεKVdt        # Calculate σ_KV explicitly at this step
        # Stress in the spring of the Kelvin-Voigt element (elastic part)
        σ_spring[i] = E2 * ε_KV[i]
    end
    return t, σ, σ_spring
end

viscous1 = LinearViscosity(5.0e19)
viscous2 = LinearViscosity(1.0e20)
elastic = Elasticity(1.0e10, 3.0e10)
elastic1 = Elasticity(1.0e10, 4.0e10)

c, x, vars, args, others = let
    # Burger's model
    #      elastic - viscous -    parallel
    #                                |
    #                   elastic --- viscous

    p = ParallelModel(viscous2, elastic)
    viscous3 = LinearViscosity(1.0e21)
    c = SeriesModel(viscous3, elastic1, p)

    vars = (; ε = 1.0e-15, θ = 1.0e-20)         # input variables (constant)
    args = (; τ = 2.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (1.0, 2.0), P0 = (0.0, 0.1))       # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

η1 = 2 * c.leafs[1].η
G1 = 2 * c.leafs[2].G
if length(c.branches) > 0
    η2 = 2 * c.branches[1].leafs[1].η
    G2 = 2 * c.branches[1].leafs[2].G
else
    η2 = 0.0
    G2 = 0.0
end

let
    function figure()
        dt = 1e9 .* [1.0, 1/2, 1/4, 1/8]
        nt = 25 .* [1.0, 2, 4, 8]
        ϵ  = zero(dt)

        for it in eachindex(dt)
            # Burgers model, numerics
            t_v, τ1, τ2, P1, P2, x1 = stress_time(c, vars, x; ntime = Int64(nt[it]), dt = dt[it]);
            t_anal, τ1_anal, τ2_anal = simulate_series_Burgers_model(G1, η1, G2, η2, vars.ε, Int64(nt[it]),  dt[it]);
            ϵ[it] = maximum(abs.(τ1 .- τ1_anal))

            # Order
            θ      = log(ϵ[1]/ϵ[2]) / log(dt[1]/dt[2])
            dt_arr = LinRange(dt[end], dt[1], 100)
            ϵ_arr  = 10 .^(log10.(ϵ[1]) .- θ.*( log10.(1 ./ (dt_arr)) .- log10.(1 ./ (dt_arr[end]))))

            SecYear = 3600 * 24 * 365.25
            fig = Figure(fontsize = 30, size = (800, 600))
            ax1 = Axis(fig[1, 1], title = L"$$Burgers model", xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]")

            lines!(ax1, t_anal / SecYear / 1.0e3, τ1_anal / 1.0e6, label = "analytical", linewidth = 5, color = :black)
            scatter!(ax1, t_v[1:10:end] / SecYear / 1.0e3, τ1[1:10:end] / 1.0e6, label = "numerical", color = :red, markersize = 15)
            #scatter!(ax1,t_v/SecYear/1e3,τ2/1e6, label="τ2")
            axislegend(ax1, position = :rb, labelsize=18)

            ax2 = Axis(fig[2, 1], title = L"$$Convergence", xlabel = L"$\log_{10}$ $\frac{1}{dt}$ [1/s]", ylabel = L"$\log_{10}$ $ϵ$ [Pa]")
            lines!(ax2, log10.(1 ./ dt_arr), log10.(ϵ_arr), color=:black, label="1st order")
            scatter!(ax2, log10.(1 ./ dt), log10.(ϵ), color=:black, label="numerics")
            axislegend(labelsize=18)

            save("docs/assets/Burgers_model.png", fig)
            display(fig)
        end
        @show ϵ
    end
    with_theme(figure, theme_latexfonts())
end

