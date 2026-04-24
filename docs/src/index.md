# RheologyCalculator.jl

`RheologyCalculator.jl` builds and solves local rheological models from small
viscous, elastic, and plastic building blocks. Rheologies can be composed in
series, in parallel, or in nested hybrid networks, then converted into a
nonlinear residual system solved with Newton iterations.

The typical workflow is:

1. Define rheological elements, such as viscosity, elasticity, or plasticity.
2. Compose them with [`SeriesModel`](@ref) and [`ParallelModel`](@ref).
3. Build an initial solver vector with [`initial_guess_x`](@ref).
4. Solve with [`solve`](@ref).
5. Update elastic history with
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

See [Composites](@ref composites) for model construction and [API](@ref) for the public
package interface.
