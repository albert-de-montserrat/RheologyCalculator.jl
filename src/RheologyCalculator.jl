"""
    RheologyCalculator

Build and solve local rheological models assembled from viscous, elastic, and
plastic elements.

The package core provides the composite containers (`SeriesModel`,
`ParallelModel`), equation generation, solver utilities, and the state-function
interface that concrete rheologies extend. Example rheology definitions live in
the repository's `rheologies/` directory and can be used as templates for
application-specific material laws.
"""
module RheologyCalculator

using DifferentiationInterface, StaticArrays, LinearAlgebra
import ForwardDiff: ForwardDiff

import Base.IteratorsMD.flatten

include("rheology_types.jl")
export AbstractViscosity, AbstractPlasticity, AbstractElasticity

include("state_functions.jl")

include("composite.jl")
export CompositeModel, SeriesModel, ParallelModel

include("kwargs.jl")

include("equations.jl")
export generate_equations

include("others.jl")

include("post_calculations.jl")

include("initial_guess.jl")
export initial_guess_x, x_keys

include("normalize_x.jl")
export normalisation_x

include("solver.jl")
export solve

include("strain_rate_correction.jl")
export effective_strain_rate_correction

include("print_rheology.jl")

end # module Rheology
