# Rheologies

Concrete rheology elements are regular Julia types that subtype one of the
abstract rheology supertypes:

```@docs
RheologyCalculator.AbstractRheology
RheologyCalculator.AbstractViscosity
RheologyCalculator.AbstractElasticity
RheologyCalculator.AbstractPlasticity
```

Each concrete element declares the state functions it contributes in series and
parallel composition. The state-function interface is documented below and in
[API](@ref).

## Element Interface

A concrete rheology is a small immutable Julia type plus methods for:

- `series_state_functions(r)`: functions used when `r` sits in a
  [`SeriesModel`](@ref).
- `parallel_state_functions(r)`: functions used when `r` sits in a
  [`ParallelModel`](@ref).
- State functions such as `compute_strain_rate`, `compute_stress`,
  `compute_pressure`, or plastic consistency functions.
- `history_kwargs(r)`, when values from `others` should be indexed per element.

For example, a deviatoric Newtonian viscosity contributes strain rate in series
and stress in parallel:

```julia
struct LinearViscosity{T} <: AbstractViscosity
    η::T
end

RheologyCalculator.series_state_functions(::LinearViscosity) =
    (RheologyCalculator.compute_strain_rate,)
RheologyCalculator.parallel_state_functions(::LinearViscosity) =
    (RheologyCalculator.compute_stress,)

RheologyCalculator.compute_strain_rate(r::LinearViscosity; τ = 0, kwargs...) = τ / (2 * r.η)
RheologyCalculator.compute_stress(r::LinearViscosity; ε = 0, kwargs...) = 2 * r.η * ε
```

The solver passes local arguments as keywords. Unknowns come from the solver
vector `x`, prescribed inputs come from `vars`, and auxiliary values come from
`others`. Element-local history fields are selected with
[`history_kwargs`](@ref RheologyCalculator.history_kwargs); for elastic elements
the default history fields are `τ0` and `P0`, and for viscous elements the
default is `d`.

## Creep Laws

### Linear Viscosity

```julia
LinearViscosity(η)
```

For scalar deviatoric quantities:

```math
\tau = 2\eta\dot\varepsilon
```

### Power-Law Viscosity

```julia
PowerLawViscosity(η, n)
```

The current implementation uses:

```math
\dot\varepsilon = \frac{\tau^n}{2\eta}
```

### Diffusion Creep

```julia
DiffusionCreep(n, r, p, A, E, V, R)
```

where `n` is the stress exponent, `r` is the water-fugacity exponent, `p` is the
grain-size exponent, `A` is the prefactor, `E` is the activation energy, `V` is
the activation volume, and `R` is the gas constant.

### Dislocation Creep

```julia
DislocationCreep(n, r, A, E, V, R)
```

where `n` is the stress exponent, `r` is the water-fugacity exponent, `A` is the
prefactor, `E` is the activation energy, `V` is the activation volume, and `R`
is the gas constant.

## Elasticity

### Compressible Elasticity

```julia
Elasticity(G, K)
```

where `G` and `K` are the shear and bulk moduli.

### Incompressible Elasticity

```julia
IncompressibleElasticity(G)
```

where `G` is the shear modulus.

### Bulk Elements

```julia
BulkElasticity(K)
BulkViscosity(χ)
```

## Plastic Failure

### Drucker-Prager

```julia
DruckerPrager(C, ϕ, ψ)
```

where `C` is the cohesion, and `ϕ` and `ψ` are the friction and dilation angles.

## State Functions

Rheologies extend these methods to participate in equation generation:

- [`compute_strain_rate`](@ref RheologyCalculator.compute_strain_rate)
- [`compute_stress`](@ref RheologyCalculator.compute_stress)
- [`compute_volumetric_strain_rate`](@ref RheologyCalculator.compute_volumetric_strain_rate)
- [`compute_pressure`](@ref RheologyCalculator.compute_pressure)
- [`compute_lambda`](@ref RheologyCalculator.compute_lambda)
- [`compute_lambda_parallel`](@ref RheologyCalculator.compute_lambda_parallel)
- [`compute_plastic_strain_rate`](@ref RheologyCalculator.compute_plastic_strain_rate)
- [`compute_volumetric_plastic_strain_rate`](@ref RheologyCalculator.compute_volumetric_plastic_strain_rate)
- [`compute_plastic_stress`](@ref RheologyCalculator.compute_plastic_stress)
- [`compute_viscosity`](@ref RheologyCalculator.compute_viscosity)
- [`compute_viscosity_series`](@ref RheologyCalculator.compute_viscosity_series)
- [`compute_viscosity_parallel`](@ref RheologyCalculator.compute_viscosity_parallel)
