using RheologyCalculator, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie
import Statistics: mean

include("../rheologies/RheologyDefinitions.jl")

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))
    
second_invariant(ε) = sqrt((ε[1]^2+ε[2]^2)/2 + ε[3]^2)

function stress_time(c, ε, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    t_v  = zeros(ntime)
    τe   = (0e0,), (0e0,), (0e0,)
    P_e  = (0.0,)
    t    = 0.0
    εII  = second_invariant(ε)
    for i in 2:ntime
        others = (; dt = dt, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
        x   = solve(c, x, ε, τe, vars, others, verbose = false)
        τII = x[1]
        ε_corr = effective_strain_rate_correction(c, ε, τe, others)
        εeff = ε .+ ε_corr
        εII_eff = second_invariant(εeff)
        η   = τII / (2 * εII_eff)

        # @show x, ε_corr, τe
        τxx, τyy, τxy   = @. 2 * εeff * η 
        # τxx_e = compute_stress_elastic(c, SA[τxx], merge(others, (; τ0 = τe[1])))
        # τyy_e = compute_stress_elastic(c, SA[τyy], merge(others, (; τ0 = τe[2])))
        # τxy_e = compute_stress_elastic(c, SA[τxy], merge(others, (; τ0 = τe[3])))
        τe    = (τxx...,), (τyy...,), (τxy...,)
    
        τ1[i] = τII
        t += others.dt
        τ_an[i] = analytical_solution(εII, t, c.leafs[2].G, c.leafs[1].η)
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
    
    ε      = 1.0e-14, -1.0e-14, 0e0
    τ0     = (0e0,), (0e0,), (0e0,)
    vars   = (; θ = 1.0e-20)  # input variables (constant)
    args   = (; τ = 2.0e6, P = 1.0e6)      # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others) .* 1e6 

    c, x, ε, τ0, vars, args, others
end

let
    function figure()
        dt = 1e10 .* [1.0, 1/2, 1/4, 1/8]
        nt = 1_000 .* [1.0, 2, 4, 8]
        ϵ  = zero(dt) 

        it = 4
        # for it in eachindex(dt)
            t_v, τ, τ_an = stress_time(c, ε, vars, x; ntime = Int64(nt[it]), dt = dt[it])
            ϵ[it] = maximum(abs.(τ .- τ_an))

            # Order
            θ      = log(ϵ[1]/ϵ[2]) / log(dt[1]/dt[2])
            dt_arr = LinRange(dt[end], dt[1], 100)
            ϵ_arr  = 10 .^(log10.(ϵ[1]) .- θ.*( log10.(1 ./ (dt_arr)) .- log10.(1 ./ (dt_arr[end]))))

            SecYear = 3600 * 24 * 365.25
            fig = Figure(fontsize = 30, size = (800, 600) .* 1)

            ax1 = Axis(fig[1, 1], title = L"$$Maxwell model", xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]")
            lines!(ax1, t_v / SecYear / 1.0e3, τ_an / 1.0e6, color=:black, label = "analytical")
            scatter!(ax1, t_v / SecYear / 1.0e3, τ / 1.0e6,  color=:red, label = "numerical")
            # scatter!(ax1, t_v[1:1000:end] / SecYear / 1.0e3, τ[1:1000:end] / 1.0e6,  color=:red, label = "numerical")
            axislegend(ax1, position = :rb, labelsize=18)
            
            ax2 = Axis(fig[2, 1], title = L"$$Convergence", xlabel = L"$\log_{10}$ $\frac{1}{dt} $ [1/s]", ylabel = L"$\log_{10}$ $ϵ$ [Pa]")
            lines!(ax2, log10.(1 ./ dt_arr), log10.(ϵ_arr), color=:black, label="1st order")
            scatter!(ax2, log10.(1 ./ dt), log10.(ϵ), color=:black, label="numerics")
            axislegend(labelsize=18)

            save("docs/assets/Maxwell_VE_model.png", fig)
            display(fig)
        # end
    end
    with_theme(figure, theme_latexfonts())
end

# (x, ε_corr, τe) = ([1.9801980198019802], (7.001057239470769e-21, -7.001057239470769e-21, 0.0), ((1.4002114478941536,), (-1.4002114478941536,), (0.0,)))

# N     = 1000
# τ_num = zeros(N)
# τ_an  = zeros(N)
# τo    = 1e3,-1e3,0e0
# εII   = second_invariant(ε)
# t = 0.0
# dt = 1e10
# for i in 2:N
#     ε_corr  = @. 1 / elastic.G/2/others.dt * τo
#     ε_eff   = ε .+ ε_corr
#     εII_eff = second_invariant(ε_eff)
#     η_eff   = 1/(1 / (viscous.η)  + 1 / elastic.G/others.dt)
#     τII     = 2 * εII_eff * η_eff

#     η       = τII / (2 * εII_eff)

#     # @show x, ε_corr, τe
#     τo = τxx, τyy, τxy   = @. 2 * ε_eff * η 
#     τ_an[i] = analytical_solution(εII, t, elastic.G, viscous.η)
#     τ_num[i] = τII

#     t += others.dt
# end
# f,ax,=lines(τ_an./1e6)
# lines!(ax,τ_num./1e6)
# f

# ε_corr  = @. 1 / elastic.G/2/others.dt * τo
# ε_eff   = ε .+ ε_corr

# a = τ_num[1]
# εII_eff = second_invariant(ε_eff)
# J = (1 / (viscous.η) /2  + 1 / elastic.G/others.dt/2)
# r = -J  * x[1] + εII_eff
# J\-r

# RC.compute_residual(c, x, vars, others)   # initial residual

# eqs      = RC.generate_equations(c)
# args_all = RC.generate_args_template(eqs, x, others)
# # evaluates the self-components of the residual
# residual1 = RC.evaluate_state_functions(eqs, args_all, others)
# residual2 = RC.add_children(residual1, x, eqs)
# residual3 = RC.subtract_parent(residual2, x, eqs, vars)
# @edit RC.subtract_parent(residual2, x, eqs, vars)
# r

# a=RC.compute_strain_rate(viscous; args_all[1]...)
# b=RC.compute_strain_rate(elastic; args_all[1]...)
# d = a+b 
# d == residual1[1]
# d - εII