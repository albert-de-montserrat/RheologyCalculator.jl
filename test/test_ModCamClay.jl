using LinearAlgebra, Statistics, StaticArrays

include("../rheologies/ModCamClay.jl")

@testset "Mod. Cam-Clay" begin
    function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
        # Extract elastic stresses/pressure from solution vector
        τ1      = zeros(ntime)
        P1      = zeros(ntime)
        t_v     = zeros(ntime)
        τ_e     = (0.0,)
        P_e     = (0.0e6,)
        P1[1]   = P_e[1]
        τ1[1]   = τ_e[1]
        x       = SA[τ1[1],0, P1[1]]
        t       = 0.0
        for i in 2:ntime
            others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
            
            x = RheologyCalculator.solve(c, x, vars, others, verbose = false, xnorm=xnorm)
            
            t += others.dt
            
            τ_e = compute_stress_elastic(c, x, others)
            P_e = compute_pressure_elastic(c, x, others)
            τ1[i] = τ_e[1]
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
        vars = (; ε = 0*7.0e-14, θ = 7.0e-15)
        # guess variables (we solve for these, differentiable)
        args = (; τ = 0.0e3, P = 0.3e6, λ = 0)
        # other non-differentiable variables needed to evaluate the state functions
        others = (; dt = 1.0e5, τ0 = (0e0, ), P0 = (0.3e6, ))

        x       = initial_guess_x(c, vars, args, others)
        char_τ  = plastic.r*100
        char_ε  = vars.ε+vars.θ
        xnorm   = normalisation_x(c, char_τ, char_ε)

        c, x, xnorm, vars, args, others
    end

    SecYear = 3600 * 24 * 365.25
    t_v, τ, P = stress_time(c, (; ε = 0*7.0e-14, θ =   7.0e-15), x, xnorm, others; ntime = 11, dt = SecYear*2)
    @test mean(τ) ≈ 0.0
    @test mean(P) ≈ -89851.02545454343
    @test any(!isnan, τ)
    @test any(!isnan, P)

    t_v, τ, P = stress_time(c, (; ε =   0*7.0e-14, θ = -7.0e-15), x, xnorm, others; ntime = 1300, dt = 1e8)
    @test mean(τ) ≈ 0.0
    @test mean(P) ≈  8.39495381358509e7
    @test any(!isnan, τ)
    @test any(!isnan, P)

    t_v, τ, P = stress_time(c, (; ε =   7.0e-14, θ =   -4.0e-15), x, xnorm, others; ntime = 700, dt = 1e8)
    @test mean(τ) ≈ 4.786269881124213e7
    @test mean(P) ≈ 3.0553535284802888e7
    @test any(!isnan, τ)
    @test any(!isnan, P)

end