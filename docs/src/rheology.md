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

```@docs
RheologyCalculator.compute_strain_rate
RheologyCalculator.compute_stress
RheologyCalculator.compute_volumetric_strain_rate
RheologyCalculator.compute_pressure
RheologyCalculator.compute_lambda
RheologyCalculator.compute_lambda_parallel
RheologyCalculator.compute_plastic_strain_rate
RheologyCalculator.compute_volumetric_plastic_strain_rate
RheologyCalculator.compute_plastic_stress
RheologyCalculator.compute_viscosity
RheologyCalculator.compute_viscosity_series
RheologyCalculator.compute_viscosity_parallel
```
