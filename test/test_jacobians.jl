using Test, RheologyCalculator, ForwardDiff, StaticArrays
import RheologyCalculator: compute_residual


viscous1    = LinearViscosity(5.0e19)
viscous2    = LinearViscosity(1.0e20)
viscousbulk = BulkViscosity(1.0e18)
powerlaw    = PowerLawViscosity(5.0e19, 3)
drucker     = DruckerPrager(1.0e6, 10.0, 0.0)
elastic     = Elasticity(1.0e10, 1.0e12)
elastic1    = Elasticity(1.0e11, 1.0e13)

elasticbulk = BulkElasticity(1.0e10)
elasticinc  = IncompressibleElasticity(1.0e10)

diffusion   = DiffusionCreep(1, 1, 1, 1.5e-3, 1, 1, 1)
dislocation = DislocationCreep(3.5, 1, 1.1e-16, 1, 1, 1)

@testset "Jacobian 1" begin
    # elastic - viscous
    c      = SeriesModel(elastic, viscous1)
    vars   = (; ε = 1.0e-15, θ = 1.0e-20) # input variables (constant)
    args   = (; τ = 1.0e2, P = 1.0e6)     # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10)              # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
    ]

    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        1.5e-20   0.0
        0.0       1.0e-22
    ]
    @test J ≈ J_sol atol = 1.0e-15
end


@testset "Jacobian 2" begin
    # elastic - viscous -- parallel
    #                         |
    #                viscous --- viscous
    s1     = SeriesModel(viscous1, viscous2)
    p      = ParallelModel(viscous1, viscous2)
    c      = SeriesModel(elastic1, viscous1, p)
    vars   = (; ε = 1.0e-15, θ = 1.0e-20)      # input variables (constant)
    args   = (; τ = 1.0e3, P = 1.0e6) # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10)          # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)[1], # local  guess(es)
    ]

    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        1.5e-20   1.0      0.0
        -1.0      3.0e20   0.0
        0.0       0.0      1.0e-22
    ]
    @test J ≈ J_sol atol = 1.0e-15
end


@testset "Jacobian 3" begin

    s1     = SeriesModel(viscous1, viscous2)
    p      = ParallelModel(s1, viscous2)
    c      = SeriesModel(viscous1, p)
    vars   = (; ε = 1.0e-15) # input variables (constant)
    args   = (; τ = 1.0e2)   # guess variables (we solve for these, differentiable)
    others = (;)             # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
        values(args)..., # local  guess(es)
    ]

    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        1.0e-20   1.0      0.0
        -1.0      2.0e20   1.0
        0.0      -1.0      1.5e-20
    ]

    @test J ≈ J_sol atol = 1.0e-15
end

@testset "Jacobian " begin
    # viscous -- parallel
    #               |
    #      viscous --- viscous
    p      = ParallelModel(viscous1, viscous2)
    c      = SeriesModel(viscous1, p)
    vars   = (; ε = 1.0e-15) # input variables (constant)
    args   = (; τ = 1.0e2)   # guess variables (we solve for these, differentiable)
    others = (;)             # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
    ]

    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        1.0e-20   1.0
        -1.0      3.0e20
    ]

    @test J ≈ J_sol atol = 1.0e-15
end

@testset "Jacobian " begin
    # viscous -- parallel
    #               |
    #      viscous --- viscous
    #         |
    #      viscous
    s1     = SeriesModel(viscous1, viscous2)
    p      = ParallelModel(s1, viscous2)
    c      = SeriesModel(viscous1, p)
    vars   = (; ε = 1.0e-15) # input variables (constant)
    args   = (; τ = 1.0e2)   # guess variables (we solve for these, differentiable)
    others = (;)             # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
        values(args)..., # local  guess(es)
    ]
    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        1.0e-20   1.0      0.0
        -1.0      2.0e20   1.0
        0.0      -1.0      1.5e-20
    ]

    @test J ≈ J_sol atol = 1.0e-15
end

@testset "Jacobian " begin
    #           parallel
    #               |
    #      viscous --- viscous
    #         |
    #      viscous
    s1     = SeriesModel(viscous1, viscous2)
    c      = ParallelModel(viscous1, s1) |> SeriesModel

    vars   = (; ε = 1.0e-15) # input variables (constant)
    args   = (; τ = 1.0e2)   # guess variables (we solve for these, differentiable)
    others = (;)             # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        0.0      1.0      0.0
        -1.0     1.0e20   1.0
        0.0     -1.0      1.5e-20
    ]

    @test J ≈ J_sol atol = 1.0e-15
end

@testset "Jacobian " begin
    # viscous -- parallel    --      parallel
    #               |                   |
    #      viscous --- viscous  viscous --- viscous
    #         |
    #      viscous
    s1     = SeriesModel(viscous1, viscous2)
    p      = ParallelModel(s1, viscous2)
    p1     = ParallelModel(viscous1, viscous2)
    c      = SeriesModel(viscous1, p, p1)
    vars   = (; ε = 1.0e-15) # input variables (constant)
    args   = (; τ = 1.0e2)   # guess variables (we solve for these, differentiable)
    others = (;)             # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)
    J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
    J_sol = SA[
        1.0e-20   1.0      1.0       0.0
        -1.0      2.0e20   1.0       0.0
        0.0      -1.0      1.5e-20   0.0
        -1.0      0.0      0.0       3.0e20
    ]

    @test J ≈ J_sol atol = 1.0e-15
end