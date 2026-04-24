using LinearAlgebra

@testset "VE Model " begin
    analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))

    function stress_time(c, vars, x, xnorm; ntime = 200, dt = 1.0e8)
        # Extract elastic stresses/pressure from solution vector
        τ1   = zeros(ntime)
        τ_an = zeros(ntime)
        t_v  = zeros(ntime)
        τ_e  = (zero_stress_tensor_2D(),)
        P_e  = (0.0,)
        t    = 0.0
        εII  = second_invariant_2D(vars.ε)
        for i in 2:ntime
            others = (; dt = dt, τ0 = τ_e, P0 = P_e)

            x      = solve(c, x, vars, others, verbose = false, xnorm0=xnorm)
            τ1[i]  = x[1]
            t     += others.dt
            τ_an[i] = analytical_solution(εII, t, c.leafs[2].G, c.leafs[1].η)
            τ_e     = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
        
            t_v[i] = t
        end

        return t_v, τ1, τ_an
    end

    c, x, xnorm, vars, args, others = let

        viscous = LinearViscosity(1e22)
        elastic = IncompressibleElasticity(10e9)

        # Maxwell viscoelastic model
        # elastic --- viscous

        c      = SeriesModel(viscous, elastic)

        vars   = vars_2D(1.0e-14, 1.0e-20)                # input variables (constant)
        args   = (; τ = 2.0e3, P = 1.0e6)                    # guess variables (we solve for these, differentiable)
        others = (; dt = 1.0e10, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

        x = initial_guess_x(c, vars, args, others)

        char_τ = 1e6
        char_ε = second_invariant_2D(vars.ε) + abs(vars.θ)
        xnorm  = normalisation_x(c, char_τ, char_ε)

        c, x, xnorm, vars, args, others
    end

    _, τ, τ_an = stress_time(c, vars, x, xnorm; ntime = 10_000, dt = 1e9)

    @test error = norm((abs.(τ_an .- τ) ./ τ_an)[2:end]) < 0.015
    @test all(isfinite, τ)
end
