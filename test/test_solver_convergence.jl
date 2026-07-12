using Test, StaticArrays
import RheologyCalculator: mynorm, _direct_leaf_elastic_correction, second_invariant_value

function converged_residual(c, x, vars0, others, xnorm)
    ε_corr = _direct_leaf_elastic_correction(c, vars0.ε, others)
    εII    = second_invariant_value(vars0.ε .+ ε_corr)
    vars   = merge(vars0, (; ε = εII))
    r = compute_residual(c, x, vars, others)
    return mynorm(r, xnorm)
end

@testset "solver convergence regression" begin
    @testset "yield-crossing hard step" begin
        viscous     = LinearViscosity(1.0e22)
        viscous_reg = LinearViscosity(1.0e20)
        elastic     = IncompressibleElasticity(10.0e9)
        plastic     = DruckerPrager(10.0e6, 30.0, 0.0)
        c = SeriesModel(viscous, elastic, ParallelModel(viscous_reg, plastic))

        vars   = vars_2D(1.0e-14)
        args   = (; τ = 2.0e3, λ = 0.0)
        others = (; dt = 1.0e8, τ0 = (zero_stress_tensor_2D(),), P = 1.0e6, P0 = (0.0,))
        x      = initial_guess_x(c, vars, args, others)
        xnorm  = normalisation_x(c, 1.0e6, second_invariant_2D(vars.ε))

        # Step the constant-strain-rate loading path forward to just before
        # the step at which the stress state crosses the yield surface.
        τ_e, P_e = (zero_stress_tensor_2D(),), (0.0,)
        x_hard, o_hard = x, others
        for _ in 1:469
            o = (; dt = others.dt, τ0 = τ_e, P = others.P, P0 = P_e)
            x_hard, o_hard = x, o
            x = solve(c, x, vars, o; xnorm0 = xnorm)
            τ_e = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, o)
        end

        x_conv = solve(c, x_hard, vars, o_hard; xnorm0 = xnorm, itermax = 10)
        @test converged_residual(c, x_conv, vars, o_hard, xnorm) < 1.0e-8
    end

    @testset "cold-start overstress" begin
        function run_case(overstress_fac, ϕ, C, εrate, dt)
            viscous     = LinearViscosity(1.0e22)
            viscous_reg = LinearViscosity(1.0e20)
            elastic     = IncompressibleElasticity(30.0e9)
            plastic     = DruckerPrager(C, ϕ, 0.0)
            c = SeriesModel(viscous, elastic, ParallelModel(plastic, viscous_reg))

            εᵢⱼ  = tensor_strain_rate_2D(εrate)
            vars = (; ε = εᵢⱼ, θ = 0.0e0)

            P            = 5.0e6
            yield_stress = P * plastic.sinϕ + plastic.C * plastic.cosϕ
            τ0ᵢⱼ         = (stress_tensor_from_invariant_2D(overstress_fac * yield_stress, εᵢⱼ),)

            args   = (; τ = 0.0, λ = 0.0, P = P)
            others = (; dt = dt, τ0 = τ0ᵢⱼ, P0 = (0.0,))

            x     = initial_guess_x(c, vars, args, others)
            xnorm = normalisation_x(c, 1.0e6, second_invariant_2D(vars.ε))

            x_conv = solve(c, x, vars, others; xnorm0 = xnorm, itermax = 10)
            return converged_residual(c, x_conv, vars, others, xnorm)
        end

        cases = (
            ("baseline",             1.1,  5.0,  10.0e6, 1.0e-14, 1.0e8),
            ("far above yield",      3.0,  5.0,  10.0e6, 1.0e-14, 1.0e8),
            ("steep friction angle", 1.1,  45.0, 10.0e6, 1.0e-14, 1.0e8),
            ("lower cohesion",       1.1,  5.0,  1.0e6,  1.0e-14, 1.0e8),
            ("small timestep",       1.1,  5.0,  10.0e6, 1.0e-14, 1.0e6),
        )

        for (label, overstress_fac, ϕ, C, εrate, dt) in cases
            @test run_case(overstress_fac, ϕ, C, εrate, dt) < 1.0e-8
        end
    end
end
