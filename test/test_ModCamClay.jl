using LinearAlgebra, Statistics, StaticArrays

include("../rheologies/ModCamClay.jl")

@testset "Mod. Cam-Clay" begin
    function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
        # Extract elastic stresses/pressure from solution vector
        τ1      = zeros(ntime)
        P1      = zeros(ntime)
        t_v     = zeros(ntime)
        τ_e     = others.τ0
        P_e     = others.P0
        P1[1]   = P_e[1]
        τ1[1]   = second_invariant_2D(τ_e[1])
        x       = SA[τ1[1], 0, P1[1]]
        t       = 0.0
        for i in 2:ntime
            others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
            
            x = RheologyCalculator.solve(c, x, vars, others, verbose = false, xnorm0=xnorm)
            
            t += others.dt
            
            τ_e = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
            P_e = compute_pressure_elastic(c, x, others)
            τ1[i] = second_invariant_2D(τ_e[1])
            P1[i] = P_e[1]

            t_v[i] = t
        end

        return t_v, τ1, P1
    end

    c, x, xnorm, vars, args, others = let

        viscous = LinearViscosity(1e23)
        elastic = Elasticity(1e10, 2e11)
        plastic = ModCamClay(; M=0.9, N=0.5, r=1e8, β=.1, Pt=-1e5, η_vp=1e20) 

        # Maxwell viscoelastic model
        c  = SeriesModel(viscous, elastic, plastic)

        # input variables (constant)
        vars = vars_2D(0*7.0e-14, 7.0e-15)
        # guess variables (we solve for these, differentiable)
        args = (; τ = 0.0e3, P = 0.3e6, λ = 0)
        # other non-differentiable variables needed to evaluate the state functions
        others = (; dt = 1.0e5, τ0 = (zero_stress_tensor_2D(),), P0 = (0.3e6, ))

        x       = initial_guess_x(c, vars, args, others)
        char_τ  = plastic.r*100
        char_ε  = second_invariant_2D(vars.ε) + abs(vars.θ)
        xnorm   = normalisation_x(c, char_τ, char_ε)

        c, x, xnorm, vars, args, others
    end

    SecYear = 3600 * 24 * 365.25
    t_v, τ, P = stress_time(c, vars_2D(0*7.0e-14, 7.0e-15), x, xnorm, others; ntime = 11, dt = SecYear*2)
    @test mean(τ) ≈ 0.0
    @test mean(P) ≈ 1489.7454545486437
    @test all(isfinite, τ)
    @test all(isfinite, P)

    t_v, τ, P = stress_time(c, vars_2D(0*7.0e-14, -7.0e-15), x, xnorm, others; ntime = 1300, dt = 1e8)
    @test mean(τ) ≈ 0.0
    @test mean(P) ≈ 8.41662039460073e7
    @test all(isfinite, τ)
    @test all(isfinite, P)

    t_v, τ, P = stress_time(c, vars_2D(7.0e-14, -4.0e-15), x, xnorm, others; ntime = 700, dt = 1e8)
    @test mean(τ) ≈ 4.7872457278776795e7
    @test mean(P) ≈ 3.080394041337929e7
    @test all(isfinite, τ)
    @test all(isfinite, P)

end
