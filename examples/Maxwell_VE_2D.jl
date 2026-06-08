using RheologyCalculator, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

# using GLMakie
import Statistics: mean

include("../rheologies/RheologyDefinitions.jl")

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))
    
second_invariant(ε) = sqrt((ε[1]^2+ε[2]^2)/2 + ε[3]^2)

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    t_v  = zeros(ntime)
    τe   = ((0e0, 0e0, 0e0),)
    P_e  = (0.0,)
    t    = 0.0
    εII  = second_invariant(vars.ε)
    ηeff = 1 / (1 / c.leafs[1].η + 1 / c.leafs[2].G / dt)
    for i in 2:ntime
        others = (; dt = dt, τ0 = τe, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
        x   = solve(c, x, vars, others, verbose = false)
        τII = x[1]
        ε_corr = effective_strain_rate_correction(c, vars.ε, τe, others)
        εeff = vars.ε .+ ε_corr
        εII_eff = second_invariant(εeff)
        η   = τII / (2 * εII_eff)

        # @show x, ε_corr, τe
        τxx, τyy, τxy   = @. 2 * εeff * η
        τxx_e = compute_stress_elastic(c, SA[τxx], merge(others, (; τ0 = τe[1][1])))
        τyy_e = compute_stress_elastic(c, SA[τyy], merge(others, (; τ0 = τe[1][2])))
        τxy_e = compute_stress_elastic(c, SA[τxy], merge(others, (; τ0 = τe[1][3])))
        τe    = ((τxx, τyy, τxy),)
    
        #τ1[i] = τII
        τ1[i] = x[1]
        
        t += others.dt
        τ_an[i] = analytical_solution(εII, t, c.leafs[2].G, c.leafs[1].η)
        t_v[i] = t
        #error("stop")
        #@show τII, εII, t   
    end

    return t_v, τ1, τ_an
end

c, x, vars, args, others = let

    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c      = SeriesModel(viscous, elastic)
    
    εᵢⱼ    = 1.0e-14, -1.0e-14, 0e0
    τ0ᵢⱼ   = ((0e0, 0e0, 0e0), )
    vars   = (; ε = εᵢⱼ, θ = 1.0e-20)  # input variables (constant)
    args   = (; τ = 2.0e6, P = 1.0e6)      # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = τ0ᵢⱼ, P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

t_v, τ, τ_an = stress_time(c, vars, x; ntime = 1000, dt = 1e10)

fig,ax,li = lines(t_v, τ_an/1e6)
scatter!(t_v, τ/1e6)

display(fig)
