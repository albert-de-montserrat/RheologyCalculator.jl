using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie, MathTeXEngine
import Statistics: mean
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    t_v  = zeros(ntime)
    τ_e  = (0.0,)
    P_e  = (0.0,)
    t   = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

        x = solve(c, x, vars, others, verbose = false)
        τ1[i] = x[1]
        t += others.dt
        τ_an[i] = analytical_solution(vars.ε, t, c.leafs[2].G, c.leafs[1].η)
        τ_e = compute_stress_elastic(c, x, others)
    
        t_v[i] = t
    end

    return t_v, τ1, τ_an
end

c, x, vars, args, others = let

    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c  = SeriesModel(viscous, elastic)

    vars   = (; ε = 1.0e-14, θ = 1.0e-20)                # input variables (constant)
    args   = (; τ = 2.0e3, P = 1.0e6)                    # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0e0, ), P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

let
    dt = 1e10 .* [1.0, 1/2, 1/4, 1/8]
    nt = 1_000 .* [1.0, 2, 4, 8]
    ϵ  = zero(dt) 

    for it in eachindex(dt)
        t_v, τ, τ_an = stress_time(c, vars, x; ntime = Int64(nt[it]), dt = dt[it])
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
        axislegend(ax, position = :rb, labelsize=18)
        
        ax2 = Axis(fig[2, 1], title = L"$$Convergence", xlabel = L"$\log_{10}$ $\frac{1}{dt} $ [1/s]", ylabel = L"$\log_{10}$ $ϵ$ [MPa]")
        lines!(ax2, log10.(1 ./ dt_arr), log10.(ϵ_arr), color=:black, label="1st order")
        scatter!(ax2, log10.(1 ./ dt), log10.(ϵ), color=:black, label="numerics")
        axislegend(labelsize=18)

        save("docs/assets/Maxwell_VE_model.png", fig)
        display(fig)
    end

end