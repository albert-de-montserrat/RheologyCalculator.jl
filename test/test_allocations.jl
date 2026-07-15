using Test

@testset "allocations" begin
    viscous1 = LinearViscosity(5.0e19)
    viscous2 = LinearViscosity(1.0e20)
    elastic  = Elasticity(1.0e10, 1.0e12)

    # flat composite: no nested branches at all
    c_flat = SeriesModel(elastic, viscous1)
    vars_flat   = (; ε = 1.0e-15, θ = 1.0e-20)
    others_flat = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    x_flat = initial_guess_x(c_flat, vars_flat, (; τ = 1.0e2, P = 1.0e6), others_flat)

    # single branch nested one level deep (pre-existing regression coverage)
    s1 = SeriesModel(viscous1, viscous2)
    p  = ParallelModel(s1, viscous2)
    c_single_branch = SeriesModel(viscous1, p)
    vars_sb   = (; ε = 1.0e-15)
    others_sb = (;)
    x_sb = initial_guess_x(c_single_branch, vars_sb, (; τ = 1.0e2), others_sb)

    # sibling branches where a non-last branch nests further: the topology
    # that exposed the equation-indexing bug in generate_offsets_parallel
    p1 = ParallelModel(viscous1, viscous2)
    c_multi_branch = SeriesModel(viscous1, p, p1)
    vars_mb   = (; ε = 1.0e-15)
    others_mb = (;)
    x_mb = initial_guess_x(c_multi_branch, vars_mb, (; τ = 1.0e2), others_mb)

    cases = (
        (c_flat, vars_flat, others_flat, x_flat),
        (c_single_branch, vars_sb, others_sb, x_sb),
        (c_multi_branch, vars_mb, others_mb, x_mb),
    )

    @testset "generate_equations" begin
        for (c, _, _, _) in cases
            generate_equations(c)
            @test (@allocated generate_equations(c)) == 0
        end
    end

    @testset "compute_residual" begin
        for (c, vars, others, x) in cases
            compute_residual(c, x, vars, others)
            @test (@allocated compute_residual(c, x, vars, others)) == 0
        end
    end

    @testset "solve" begin
        for (c, vars, others, x) in cases
            solve(c, x, vars, others)
            @test (@allocated solve(c, x, vars, others)) == 0
        end
    end
end
