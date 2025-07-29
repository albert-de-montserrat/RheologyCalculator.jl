# Rheologies

## Creep laws
### Linear viscosity

$\tau = 2\eta\dot\varepsilon$

```julia
LinearViscosity(η)
```

### Power law viscosity

$\tau = 2\eta\dot\varepsilon^n$

```julia
PowerLawViscosity(η, n)
```

### Diffusion Creep

$\dot{\varepsilon} = A \tau d^{\mathrm{p}} f_\mathrm{H_2O}^r \exp\left(-\frac{E+PV}{RT}\right)$

```julia
DiffusionCreep(n, r, A, E, V, R)
```

where 'n' is the power-law exponent,'r' is the exponent of water-fugacity,'A' is the material specific rheological parameter,'E' is the activation energy,'V' is the activation volume,'R' is the universal gas constant.

### Dislocation Creep

$\dot{\varepsilon} = A \tau d^{\mathrm{p}} f_\mathrm{H_2O}^r \exp\left(-\frac{E+PV}{RT}\right)$

```julia
DislocationCreep(n, r, A, E, V, R)
```

where 'n' is the power-law exponent,'r' is the exponent of water-fugacity,'A' is the material specific rheological parameter,'E' is the activation energy,'V' is the activation volume,'R' is the universal gas constant.

## Elasticity

### Compressible elasticity

```julia
Elasticity(G, K)
```

where 'G' and 'K' are the shear and bulk modulus, respectively.

### Incompressible Elasticity

```julia
IncompressibleElasticity(G)
```

where 'G' is the shear modulus.

## Plastic failure

### Drucker-Prager

```julia
DruckerPrager(C, ϕ, ψ)
```

where 'C' is the cohesion, and `ϕ` and `ψ` are the friction and dilation angles, respectively.


# Others
"""@docs
BulkViscosity(η)
BulkElasticity(K)
""