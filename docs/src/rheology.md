# Rheology reference

This page documents all built-in rheological elements provided in
`rheologies/RheologyDefinitions.jl`.  Each element must implement:

1. A Julia `struct` subtyping one of `AbstractViscosity`, `AbstractElasticity`, or
   `AbstractPlasticity`.
2. The `series_state_functions` and (optionally) `parallel_state_functions` methods,
   which tell the equation-generator which state functions govern the element.
3. One or more `compute_*` state functions that evaluate the constitutive relationship.

---

## Viscous elements

### `LinearViscosity`

Newton's law of viscosity:

$$\dot\varepsilon = \frac{\tau}{2\eta}, \qquad \tau = 2\eta\dot\varepsilon$$

```julia
LinearViscosity(η)
```

| Parameter | Description |
|:---|:---|
| `η` | Dynamic shear viscosity (Pa·s) |

**Example**
```julia
viscous = LinearViscosity(1e22)   # η = 1e22 Pa·s
```

---

### `PowerLawViscosity`

Power-law viscosity:

$$\dot\varepsilon = \frac{\tau^n}{2\eta}, \qquad \tau = (2\eta\,\dot\varepsilon)^{1/n}$$

```julia
PowerLawViscosity(η, n)
```

| Parameter | Description |
|:---|:---|
| `η` | Consistency parameter |
| `n` | Power-law exponent (`n=1` recovers `LinearViscosity`) |

**Example**
```julia
rock = PowerLawViscosity(1e15, 3)
```

---

### `LTPViscosity`

Low-temperature plasticity (Peierls creep):

$$\dot\varepsilon = \varepsilon_0\,\sinh\!\left(\frac{Q(\tau - \sigma_b)}{\sigma_r}\right)$$

```julia
LTPViscosity(ε0, Q, σb, σr)
```

| Parameter | Description |
|:---|:---|
| `ε0` | Reference strain rate (s⁻¹) |
| `Q`  | Activation barrier (dimensionless) |
| `σb` | Peierls/brittle strength (Pa) |
| `σr` | Reference stress (Pa) |

---

### `DiffusionCreep`

$$\dot\varepsilon = A\,\tau^n\,d^{-p}\,f_{\mathrm{H_2O}}^r \exp\!\left(-\frac{E+PV}{RT}\right)$$

```julia
DiffusionCreep(n, r, p, A, E, V, R)
```

| Parameter | Description |
|:---|:---|
| `n` | Stress exponent |
| `r` | Water-fugacity exponent |
| `p` | Grain-size exponent |
| `A` | Pre-exponential factor |
| `E` | Activation energy (J/mol) |
| `V` | Activation volume (m³/mol) |
| `R` | Universal gas constant (J/(mol·K)) |

---

### `DislocationCreep`

$$\dot\varepsilon = A\,\tau^n\,f_{\mathrm{H_2O}}^r \exp\!\left(-\frac{E+PV}{RT}\right)$$

```julia
DislocationCreep(n, r, A, E, V, R)
```

| Parameter | Description |
|:---|:---|
| `n` | Power-law exponent |
| `r` | Water-fugacity exponent |
| `A` | Pre-exponential factor |
| `E` | Activation energy (J/mol) |
| `V` | Activation volume (m³/mol) |
| `R` | Universal gas constant (J/(mol·K)) |

---

### `BulkViscosity`

$$\dot\theta = P/\chi, \qquad P = \chi\,\dot\theta$$

```julia
BulkViscosity(χ)
```

| Parameter | Description |
|:---|:---|
| `χ` | Bulk viscosity (Pa·s) |

---

## Elastic elements

Elastic elements require history variables — the stress and pressure at the end of the
**previous** time step — provided through `others` as `τ0` (deviatoric tuple) and `P0`
(volumetric tuple).

### `IncompressibleElasticity`

Purely deviatoric (incompressible) Maxwell spring:

$$\dot\varepsilon^e = \frac{\tau - \tau_0}{2G\,\Delta t}, \qquad \tau = \tau_0 + 2G\,\Delta t\,\dot\varepsilon$$

```julia
IncompressibleElasticity(G)
```

| Parameter | Description |
|:---|:---|
| `G` | Shear modulus (Pa) |

**Example**
```julia
elastic = IncompressibleElasticity(10e9)
others  = (; dt = 1e10, τ0 = (0.0,), P0 = (0.0,))
```

---

### `Elasticity`

Compressible element with deviatoric and volumetric contributions:

$$\tau = \tau_0 + 2G\,\Delta t\,\dot\varepsilon^e, \qquad P = P_0 - K\,\Delta t\,\dot\theta^e$$

```julia
Elasticity(G, K)
```

| Parameter | Description |
|:---|:---|
| `G` | Shear modulus (Pa) |
| `K` | Bulk modulus (Pa) |

Use `Elasticity` (rather than `IncompressibleElasticity`) whenever the dilatancy angle
`ψ ≠ 0` in the associated plastic element.

---

### `BulkElasticity`

Purely volumetric elastic element:

$$P = P_0 - K\,\Delta t\,\dot\theta^e$$

```julia
BulkElasticity(K)
```

| Parameter | Description |
|:---|:---|
| `K` | Bulk modulus (Pa) |

---

## Plastic elements

### `DruckerPrager`

Pressure-dependent Drucker-Prager yield criterion:

$$F(\tau, P) = \tau - P\sin\phi - C\cos\phi \leq 0$$

Non-associative flow rule (`ψ ≠ ϕ`):

$$\dot\varepsilon^{\mathrm{pl}} = \lambda\,\partial_\tau Q, \qquad Q = \tau - P\sin\psi$$

```julia
DruckerPrager(C, ϕ, ψ)
```

| Parameter | Description |
|:---|:---|
| `C` | Cohesion (Pa) |
| `ϕ` | Internal friction angle (degrees) |
| `ψ` | Dilatancy angle (degrees); `ψ=0` → non-dilatant, no volumetric plastic flow |

When `ψ ≠ 0`, pair with `Elasticity(G, K)` (not `IncompressibleElasticity`).

**Non-dilatant VEP**
```julia
plastic = DruckerPrager(10e6, 30, 0)
c       = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(10e9), plastic)
```

**Dilatant VEP** (see also `examples/Maxwell_VEP_dilatant.jl`)
```julia
plastic = DruckerPrager(10e6, 30, 10)
c       = SeriesModel(LinearViscosity(1e22), Elasticity(10e9, 20e9), plastic)
```

---

## Defining a custom rheology

Adding a new creep law requires three steps:

```julia
using RheologyCalculator
import RheologyCalculator: series_state_functions, parallel_state_functions,
                            compute_strain_rate, compute_stress

# 1. Define the struct
struct MyViscosity{T} <: AbstractViscosity
    η::T
end

# 2. Register with the equation generator
@inline series_state_functions(::MyViscosity)   = (compute_strain_rate,)
@inline parallel_state_functions(::MyViscosity) = (compute_stress,)

# 3. Implement the constitutive relation
@inline compute_strain_rate(r::MyViscosity; τ = 0, kwargs...) = τ / (2 * r.η)
@inline compute_stress(r::MyViscosity; ε = 0, kwargs...)      = 2 * r.η * ε
```

After this, `MyViscosity` can be composed inside any `SeriesModel` or `ParallelModel`
without any further changes.
