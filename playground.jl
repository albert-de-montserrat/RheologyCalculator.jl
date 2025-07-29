using ForwardDiff, StaticArrays, LinearAlgebra
using GLMakie

import Base.IteratorsMD.flatten

include("rheology_types.jl")
include("state_functions.jl")
include("kwargs.jl")
include("composite.jl")
include("recursion.jl")
include("equations.jl")
include("others.jl")
include("function_utils.jl")
include("../src/print_rheology.jl")


function bt_line_search(Δx, J, x, r, composite, vars, others; α = 1.0, ρ = 0.5, c = 1.0e-4, α_min = 1.0e-8) where {N}

    perturbed_x = @. x - α * Δx
    perturbed_r = compute_residual(composite, perturbed_x, vars, others)

    while sqrt(sum(perturbed_r .^ 2)) > sqrt(sum((r + (c * α * (J * Δx))) .^ 2))
        α *= ρ
        if α < α_min
            α = α_min
            break
        end
        perturbed_x = @. x - α * Δx
        perturbed_r = compute_residual(composite, perturbed_x, vars, others)
    end
    return α
end

function solve(c, x, vars, others)
    tol = 1.0e-9
    itermax = 1.0e3
    it = 0
    er = Inf
    # Δx      = similar(x)
    local α
    while er > tol
        it += 1
        r = compute_residual(c, x, vars, others)
        J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
        Δx = J \ r
        α = bt_line_search(Δx, J, x, r, c, vars, others)
        x -= α .* Δx

        er = norm(iszero(xᵢ) ? 0.0e0 : Δxᵢ / abs(xᵢ) for (Δxᵢ, xᵢ) in zip(Δx, x)) # norm(r)

        it > itermax && break
    end
    println("Iteration: $it, Error: $er, α = $α")
    return x
end


viscous1 = LinearViscosity(5.0e19)
viscous2 = LinearViscosity(1.0e20)
viscousbulk = BulkViscosity(1.0e18)
powerlaw = PowerLawViscosity(5.0e19, 3)
drucker = DruckerPrager(1.0e6, 10.0, 0.0)
elastic = Elasticity(1.0e10, 1.0e12)
elasticbulk = BulkElasticity(1.0e10)
elasticinc = IncompressibleElasticity(1.0e10)

LTP = LTPViscosity(6.2e-13, 76, 1.8e9, 3.4e9)
diffusion = DiffusionCreep(1, 1, 1, 1.5e-3, 1, 1, 1)
dislocation = DislocationCreep(3.5, 1, 1.1e-16, 1, 1, 1)

c, x, vars, args, others = let
    # elastic - viscous -- parallel
    #                         |
    #                viscous --- viscous
    s1 = SeriesModel(viscous1, viscous2)
    p = ParallelModel(viscous1, viscous2)
    c = SeriesModel(elastic, viscous1, p)
    vars = (; ε = 1.0e-15, θ = 1.0e-20)      # input variables (constant)
    args = (; τ = 1.0e3, P = 1.0e6) # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)[1], # local  guess(es)
    ]

    c, x, vars, args, others
end

eqs = generate_equations(c)
r = compute_residual(c, x, vars, others)
J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)

ε = exp10.(LinRange(log10(1.0e-15), log10(1.0e-11), 1000))
τ = similar(ε)
# args   = (; τ = 1e10)   # guess variables (we solve for these, differentiable)
# others = (; dt = 1e10) # other non-differentiable variables needed to evaluate the state functions
for i in eachindex(ε)
    # vars = (; ε = ε[i]) # input variables (constant)
    vars = (; ε = ε[i], θ = 1.0e-15) # input variables (constant)
    sol = solve(c, x, vars, others)
    τ[i] = sol[1]
    # @show sol
end

f, ax, h = scatterlines(log10.(ε), τ)
ax.xlabel = "εII"
ax.ylabel = "τII"
f

c, x, vars, args, others = let
    # elasticbulk - viscousbulk -- parallel
    #                                 |
    #                       viscous1 --- elasticinc
    p = ParallelModel(viscous1, elastic)
    c = SeriesModel(elasticbulk, viscousbulk, p)
    vars = (; θ = 1.0e-20, ε = 1.0e-15)      # input variables (constant)
    args = (; P = 1.0e6, τ = 1.0e3) # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(vars)[1], # global guess(es), solving for these
        values(args)[1], # local  guess(es)
        values(vars)[2], # global guess(es), solving for these
        values(args)[2], # local  guess(es)
    ]
    c, x, vars, args, others
end

eqs = generate_equations(c)
r = compute_residual(c, x, vars, others)
J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)

sol = solve(c, x, vars, others)

#=
c, x, vars, args, others = let
    # viscous1  -- drucker
    c      = SeriesModel(viscous1, drucker)
    vars   = (; ε = 1e-15)      # input variables (constant)
    args   = (; τ = 1e3, λ = 0) # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)[1], # global guess(es), solving for these
        values(args)[2], # global guess(es), solving for these
    ]
    c, x, vars, args, others
end


c, x, vars, args, others = let
    # viscous1  -- drucker
    p      = ParallelModel(drucker, viscous1)    
    c      = SeriesModel(viscous1, elastic, p)
    vars   = (; θ = 1e-20, ε = 1e-15)      # input variables (constant)
    args   = (; P = 1e6, τ = 1e3, λ = 0) # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)[2], # global guess(es), solving for these
        values(args)[1], # global guess(es), solving for these
        values(args)[3], # global guess(es), solving for these
    ]
    c, x, vars, args, others
end
=#
eqs = generate_equations(c);
length(eqs)


eqs = generate_equations(c);
eqs[1].child
eqs[2].child

eqs[1].parent
eqs[2].parent

r = compute_residual(c, x, vars, others)
#J   = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
