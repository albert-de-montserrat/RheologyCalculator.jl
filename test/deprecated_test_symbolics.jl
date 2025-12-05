using Symbolics
import RheologyCalculator: compute_residual

@variables η1, η2, G1, G2, K1, K2

viscous1 = LinearViscosity(η1)
viscous2 = LinearViscosity(η2)
elastic1 = Elasticity(G1, K1)
elastic2 = Elasticity(G2, K2)
p2 = ParallelModel(viscous2, elastic2)
p1 = ParallelModel(viscous1, elastic1)

@variables ε θ τ P dt τ0 P0

@testset "Symbolic Jacobian 1" begin

    c = SeriesModel(viscous1, p1)
    vars   = (; ε = ε, θ = θ)               # input variables (constant)
    args   = (; τ = τ, P = P)               # guess variables (we solve for these, differentiable)
    others = (; dt = dt, τ0 = τ0, P0 = P0)  # other non-differentiable variables 
    x      = initial_guess_x(c, vars, args, others)

    # Compute analytical Jacobian
    J = ForwardDiff.jacobian(x -> compute_residual(c, x, vars, others), x)
    J_sym = SA[
         1 / (2η1) 1              0      0
        -1         2η1 + 2G1*dt   0      0
         0         0              0      1
         0         0             -1     -K1*dt
    ]
    @test J === J_sym
end

@testset "Symbolic Jacobian 2" begin

    c = SeriesModel(p1, p2)
    vars   = (; ε = ε, θ = θ)               # input variables (constant)
    args   = (; τ = τ, P = P)               # guess variables (we solve for these, differentiable)
    others = (; dt = dt, τ0 = τ0, P0 = P0)  # other non-differentiable variables 
    x      = initial_guess_x(c, vars, args, others)

    # Compute analytical Jacobian
    J = ForwardDiff.jacobian(x -> compute_residual(c, x, vars, others), x)
    J_sym = SA[
        0             1             1  0      0      0
       -1  2η1 + 2G1*dt             0  0      0      0
       -1             0  2η2 + 2G2*dt  0      0      0
        0             0             0  1      1      0
        0             0            -1  0 -K1*dt      0
        0             0            -1  0      0 -K2*dt
    ]
    @test J === J_sym
end

@testset "Symbolic Jacobian 3" begin

    c = SeriesModel(viscous1, elastic1, p2)
    vars   = (; ε = ε, θ = θ)               # input variables (constant)
    args   = (; τ = τ, P = P)               # guess variables (we solve for these, differentiable)
    others = (; dt = dt, τ0 = τ0, P0 = P0)  # other non-differentiable variables 
    x      = initial_guess_x(c, vars, args, others)

    # Compute analytical Jacobian
    J = ForwardDiff.jacobian(x -> compute_residual(c, x, vars, others), x)
    J_sym = SA[
    1 / (2η1) + 1 / (2G1*dt)             1             0        0
                          -1  2η2 + 2G2*dt             0        0
                           0             0   1 / (-K1*dt)       1
                           0             0            -1   -K2*dt
    ]
    @test J === J_sym
end