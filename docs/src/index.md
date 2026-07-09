# RheologyCalculator.jl

```@docs
RheologyCalculator
```

`RheologyCalculator.jl` builds and solves local rheological models from small
viscous, elastic, and plastic building blocks. Rheologies can be composed in
series, in parallel, or in nested hybrid networks, then converted into a
nonlinear residual system solved with Newton iterations.

The computational core is intentionally separate from any particular material
catalogue. The repository's `rheologies/` directory contains example concrete
laws such as `LinearViscosity`, `Elasticity`, and `DruckerPrager`; project code
can define its own types by extending the same state-function interface.

## Installation

`RheologyCalculator.jl` is registered in the Julia General registry:

```julia
using Pkg
Pkg.add("RheologyCalculator")
```

## Loading rheological elements

The package exports the composition and solver machinery ([`SeriesModel`](@ref),
[`ParallelModel`](@ref), [`solve`](@ref), [`initial_guess_x`](@ref), …) but not
the concrete constitutive elements used in the examples (`LinearViscosity`,
`Elasticity`, `DruckerPrager`, and the rest). Those are defined in
`rheologies/RheologyDefinitions.jl` in the repository, which extends the
state-function interface documented in [Rheologies](@ref), and must be loaded
explicitly:

```julia
using RheologyCalculator
include("rheologies/RheologyDefinitions.jl")
```

The `include` path is relative to the current working directory, so the examples
here assume a clone of the repository. When the package is installed with
`Pkg.add`, `include` the file by an absolute path, or copy it into your own
project.

The typical workflow is:

1. Define rheological elements, such as viscosity, elasticity, or plasticity.
2. Compose them with [`SeriesModel`](@ref) and [`ParallelModel`](@ref).
3. Provide prescribed inputs in `vars`, initial unknown estimates in `args`,
   and auxiliary/history values in `others`.
4. Build an initial solver vector with [`initial_guess_x`](@ref).
5. Solve with [`solve`](@ref).
6. Update elastic history with
   [`compute_stress_elastic`](@ref RheologyCalculator.compute_stress_elastic)
   and [`compute_pressure_elastic`](@ref RheologyCalculator.compute_pressure_elastic),
   when needed.

```julia
using RheologyCalculator

include("rheologies/RheologyDefinitions.jl")

viscous = LinearViscosity(1e22)
elastic = IncompressibleElasticity(1e10)
c = SeriesModel(viscous, elastic)

vars = (; ε = 1.0e-14, θ = 0.0)
args = (; τ = 1.0e3, P = 0.0)
others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))

x = initial_guess_x(c, vars, args, others)
x = solve(c, x, vars, others)
```

Here `vars` contains prescribed rates (`ε`, `θ`), `args` seeds the solver
unknowns (`τ`, `P`, and any branch-local unknowns), and `others` carries values
that are not differentiated by the local Newton solve (`dt`, elastic history,
grain size, temperature, pressure-dependent parameters, and similar fields).

See [Composites](@ref composites) for model construction, [Rheologies](@ref) for
the element interface, and [API](@ref) for the generated reference.
