using LinearAlgebra

@testset "VEVP Model " begin
    function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
        # Extract elastic stresses/pressure from solutio vector
        τ1  = zeros(ntime)
        t_v = zeros(ntime)
        τ_e = (zero_stress_tensor_2D(),)
        P_e = (0.0,)
        t   = 0.0
        τy  = c.branches[1].leafs[2].C * cosd(c.branches[1].leafs[2].ϕ) + others.P * sind(c.branches[1].leafs[2].ϕ)
        x   = SA[x[1], x[2], max(x[3], second_invariant_2D(vars.ε)), max(x[4], τy * 1.01)]
        for i in 2:ntime
            others   = (; dt = dt, τ0 = τ_e, P=others.P, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
            x        = solve(c, x, vars, others, verbose = false, xnorm0=xnorm)
            τ1[i]    = x[1]
            t       += others.dt
            τ_e      = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
            t_v[i]   = t
        end
        return t_v, τ1
    end

    c, x, xnorm, vars, args, others = let
        viscous     = LinearViscosity(1e22)
        viscous_reg = LinearViscosity(1e20)
        elastic     = IncompressibleElasticity(10e9)
        plastic     = DruckerPrager(10e6, 30, 0)

        # Maxwell visco-elasto-(visco-plastic) model
        p  = ParallelModel(viscous_reg, plastic)
        c  = SeriesModel(viscous, elastic, p)

        # input variables (constant)
        vars   = vars_2D(1.0e-14)
        # guess variables (we solve for these, differentiable)
        args   = (; τ = 2.0e3, λ = 0)
        # other non-differentiable variables needed to evaluate the state functions
        others = (; dt = 1.0e8, τ0 = (zero_stress_tensor_2D(),), P = 1.0e6, P0 = (0.0, ))

        x = initial_guess_x(c, vars, args, others)
        
        char_τ = 1e6
        char_ε = second_invariant_2D(vars.ε)
        xnorm  = normalisation_x(c, char_τ, char_ε)

        c, x, xnorm, vars, args, others
    end

    _, τ = stress_time(c, vars, x, xnorm, others; ntime = 1_500, dt = 1e8)
    @test maximum(τ) ≈ 1.1049753302374864e7
    @test all(isfinite, τ)
end
