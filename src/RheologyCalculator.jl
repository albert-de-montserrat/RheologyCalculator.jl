module RheologyCalculator

using ForwardDiff, StaticArrays, LinearAlgebra

import Base.IteratorsMD.flatten

include("rheology_types.jl")
#export LinearViscosity, Elasticity, LinearViscosityStress, PowerLawViscosity, DiffusionCreep, DislocationCreep, LTPViscosity, IncompressibleElasticity, BulkViscosity, BulkElasticity, DruckerPrager
export AbstractViscosity, AbstractPlasticity, AbstractElasticity

include("state_functions.jl")

include("composite.jl")
export CompositeModel, SeriesModel, ParallelModel

include("kwargs.jl")

# include("recursion.jl")

include("equations.jl")
export generate_equations

# include("residuals.jl")

include("others.jl")
export isvolumetric

include("post_calculations.jl")

include("initial_guess.jl")
export initial_guess_x

include("solver.jl")
export solve

include("print_rheology.jl")

end # module Rheology
