using RheologyCalculator
import RheologyCalculator: compute_viscosity, compute_stress_elastic, compute_pressure_elastic, effective_strain_rate_correction
import RheologyCalculator as RC

include("rheologies/RheologyDefinitions.jl")

viscous1    = LinearViscosity(1e0)
viscous_reg = LinearViscosity(1e-2)
elastic1    = Elasticity(1e0, 1e2)
plastic1    = DruckerPrager(1.6, 30, 5.0)
p           = ParallelModel(viscous_reg, elastic1)
c           = SeriesModel(viscous1, elastic1, p)
c           = SeriesModel(viscous1, elastic1)

ε      = 1e-14, -1e-14, 0e0
τ0     = ((0e0, ), (0e0, ), (0e0,))
τ0     = ((1e0, ), (-1e0, ), (0e0,))
  
# input variables (constant)
vars   = (; ε = 1.0e-14, θ = 1.0e-20)
# guess variables (we solve for these, differentiable)
args   = (; τ = 0e0, P = 1.0e6, λ = 0)
# other non-differentiable variables needed to evaluate the state functions
others = (; dt = 1.0e8, τ0 = 0e0, P0 = 0.0)
x       = initial_guess_x(c, vars, args, others)

p           = ParallelModel(viscous_reg, elastic1)
I = 1
r = c.leafs[1]

effective_strain_rate_correction(c, ε, τ0, others)

# RC.compute_viscosity_parallel(c.branches[1].leafs, ε, others)

# import RheologyCalculator: AbstractRheology, AbstractCompositeModel, iselastic

# effective_strain_rate_correction(c::AbstractRheology, ε, τ0, others, I) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others, I)

# effective_strain_rate_correction(::Val{false}, c::AbstractRheology, ε, τ0, others, I) = 0

# function effective_strain_rate_correction(::Val{true}, c::AbstractRheology, ε, τ0, others, I)
#     η = compute_viscosity(c, merge((; ε), others))
#     correction = τ0 / (2 * η)
#     return correction
# end

# effective_strain_rate_correction(RC.iselastic(c.leafs[2]), r, ε[1], τ0[1], others, I)
# effective_strain_rate_correction(r, ε[1], τ0[1], others, I)
# effective_strain_rate_correction(c, ε, τ, others)


# update_correction_index(c::AbstractRheology, I) = update_correction_index(iselastic(c), I)

# update_correction_index(::Val{false}, I) = I
# update_correction_index(::Val{true}, I)  = I + 1