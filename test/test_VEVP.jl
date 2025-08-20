using LinearAlgebra

@testset "VEVP Model " begin
    function stress_time(c, vars, x, others; ntime = 200, dt = 1.0e8)
        # Extract elastic stresses/pressure from solutio vector
        τ1  = zeros(ntime)
        t_v = zeros(ntime)
        τ_e = (0.0,)
        P_e = (0.0,)
        t   = 0.0
        for i in 2:ntime
            others   = (; dt = dt, τ0 = τ_e, P=others.P, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
            x        = solve(c, x, vars, others, verbose = false)
            τ1[i]    = x[1]
            t       += others.dt
            τ_e      = x[1] 
            t_v[i]   = t
        end
        return t_v, τ1
    end

    c, x, vars, args, others = let
        viscous     = LinearViscosity(1e22)
        viscous_reg = LinearViscosity(1e20)
        elastic     = IncompressibleElasticity(10e9)
        plastic     = DruckerPrager(10e6, 30, 0)

        # Maxwell visco-elasto-(visco-plastic) model
        p  = ParallelModel(viscous_reg, plastic)
        c  = SeriesModel(viscous, elastic, p)

        # input variables (constant)
        vars   = (; ε = 1.0e-14)
        # guess variables (we solve for these, differentiable)
        args   = (; τ = 2.0e3, λ = 0)
        # other non-differentiable variables needed to evaluate the state functions
        others = (; dt = 1.0e8, τ0 = (0e0, ), P = 1.0e6, P0 = (0.0, ))

        x = initial_guess_x(c, vars, args, others)

        c, x, vars, args, others
    end

    _, τ = stress_time(c, vars, x, others; ntime = 1_500, dt = 1e8)
    @test maximum(τ) == 1.1049696158037923e7
    @test any(!isnan, τ)
end