using LinearAlgebra

@testset "VEP Model " begin
    function analytical_solution(ϵ, t, G, η, c, ϕ, P)
        τ  = 2 * ϵ * η * (1 - exp(-G * t / η))
        τy = c*cosd(ϕ) + P*sind(ϕ)
        return τy < τ ? τy : τ
    end

    function stress_time(c, vars, x, others; ntime = 200, dt = 1.0e8)
        # Extract elastic stresses/pressure from solutio vector
        τ1   = zeros(ntime)
        τ_an = zeros(ntime)
        t_v  = zeros(ntime)
        τ_e  = 0.0
        P_e  = 0.0
        t    = 0.0
        for i in 2:ntime
            others = (; dt = dt, τ0 = τ_e, P = others.P, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

            x       = solve(c, x, vars, others, verbose = false)
            τ1[i]   = x[1]
            t      += others.dt
            τ_an[i] = analytical_solution(vars.ε, t, c.leafs[2].G, c.leafs[1].η, c.leafs[3].C, c.leafs[3].ϕ, others.P)
            τ_e     = compute_stress_elastic(c, x, others)
            t_v[i]  = t
        end

        return t_v, τ1, τ_an
    end

    c, x, vars, args, others = let

        viscous = LinearViscosity(1e22)
        elastic = IncompressibleElasticity(10e9)
        plastic = DruckerPrager(10e6, 30, 0)

        # Maxwell visco-elasto-plastic model
        # elastic --- viscous --- plastic

        c = SeriesModel(viscous, elastic, plastic)
        
        # input variables (constant)
        vars   = (; ε = 1.0e-14, θ = 1.0e-20)
        # guess variables (we solve for these, differentiable)
        args   = (; τ = 0e0, λ = 0)
        # other non-differentiable variables needed to evaluate the state functions
        others = (; dt = 1.0e8, P = 1.0e6, τ0 = 0e0, P0 = 0.0)

        x = initial_guess_x(c, vars, args, others)

        c, x, vars, args, others
    end

    _, τ, τ_an = stress_time(c, vars, x, others; ntime = 1_500, dt = 1e8)
    @test error = norm((abs.(τ_an .- τ) ./ τ_an)[2:end]) < 0.0015
    @test any(!isnan, τ)
end