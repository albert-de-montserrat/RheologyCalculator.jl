# This implements the mode1/mode2 plasticity model proposed in:
# Popov, A. A., Berlie, N., and Kaus, B. J. P.: A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: 
#   stable two-field mixed formulation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-2469, 2025.
using Test, LinearAlgebra, Statistics
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using StaticArrays

include("../rheologies/DruckerPragerCap.jl")

@testset "VEPCap Model " begin
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
        plastic = DruckerPragerCap(; C=1e6, ϕ=30.0, ψ=10.0, η_vp=0.0, Pt=-5e5) 

        c  = SeriesModel(viscous, elastic, plastic)
        
        # input variables (constant)
        vars = (; ε = 0*7.0e-14, θ = 7.0e-15)
        # guess variables (we solve for these, differentiable)
        args = (; τ = 0.0e3, P = 0.3e6, λ = 0)
        # other non-differentiable variables needed to evaluate the state functions
        others = (; dt = 1.0e5, τ0 = (0e0, ), P0 = (0.3e6, ))

        x       = initial_guess_x(c, vars, args, others)
        char_τ  = plastic.C
        char_ε  = vars.ε+vars.θ
        xnorm   = normalisation_x(c, char_τ, char_ε)

        c, x, xnorm, vars, args, others
    end


    SecYear = 3600 * 24 * 365.25
    t_v1, τ1, P1    = stress_time(c, (; ε = 0*7.0e-14, θ =   7.0e-15), x, xnorm, others; ntime = 11, dt = SecYear*2)

    @test mean(τ1) ≈ 0.0
    @test mean(P1) ≈ -347765.3818181818
    @test any(!isnan, τ1)
    @test any(!isnan, P1)

    t_v2, τ2, P2    = stress_time(c, (; ε =   7.0e-14, θ = 0*7.0e-15), x, xnorm, others; ntime = 80, dt = 1e7)
    @test mean(τ2) ≈ 543243.3665562181
    @test mean(P2) ≈  34994.05512919912
    @test any(!isnan, τ2)
    @test any(!isnan, P2)

    t_v3, τ3, P3    = stress_time(c, (; ε =   7.0e-14, θ =   7.0e-15), x, xnorm, others; ntime = 30, dt = 2e7)
    @test mean(τ3) ≈ 397000.6363506933
    @test mean(P3) ≈ -242529.90553467823
    @test any(!isnan, τ3)
    @test any(!isnan, P3)

end