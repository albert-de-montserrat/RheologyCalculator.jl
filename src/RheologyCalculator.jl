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

using StaticArrays, LinearAlgebra
import ForwardDiff: ForwardDiff

import Base.IteratorsMD.flatten

include("core/rheology_types.jl")
export AbstractViscosity, AbstractPlasticity, AbstractElasticity

include("core/state_functions.jl")

include("core/composite.jl")
export CompositeModel, SeriesModel, ParallelModel

include("core/kwargs.jl")

include("equation_system/equations.jl")
export generate_equations

include("core/others.jl")

include("postprocessing/post_calculations.jl")

include("equation_system/initial_guess.jl")
export initial_guess_x, x_keys

include("equation_system/normalize_x.jl")
export normalisation_x

include("equation_system/solver.jl")
export solve

include("postprocessing/strain_rate_correction.jl")
export effective_strain_rate_correction

include("display/print_rheology.jl")

end # module Rheology
