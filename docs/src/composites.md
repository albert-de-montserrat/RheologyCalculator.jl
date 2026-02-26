# Composite rheology models

`RheologyCalculator.jl` builds composite rheologies by assembling individual elements in
[series](https://en.wikipedia.org/wiki/Maxwell_model),
[parallel](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Voigt_material), or arbitrarily
nested hybrid configurations. The equation system and solver are generated automatically
from the topology of the composite model at compile time.

## Core types

| Type | Description |
|:---|:---|
| `SeriesModel(elements...)` | Elements connected in series (strain rates add, stresses are equal) |
| `ParallelModel(elements...)` | Elements connected in parallel (stresses add, strain rates are equal) |

Both types accept any mix of leaf rheologies (`AbstractViscosity`, `AbstractElasticity`,
`AbstractPlasticity`) and nested `SeriesModel` / `ParallelModel` objects.

## Example: Maxwell viscoelastic model (series)

A viscous damper and an elastic spring connected in series form the classic
[Maxwell model](https://en.wikipedia.org/wiki/Maxwell_model):

```
---[viscous]---[elastic]---
```

```julia
viscous = LinearViscosity(1e22)          # η = 1e22 Pa·s
elastic = IncompressibleElasticity(10e9)  # G = 10 GPa

maxwell = SeriesModel(viscous, elastic)
```

In a series model the deviatoric *stress* is the same across all elements, and the total
*strain rate* is the sum of contributions from each element.

## Example: Kelvin–Voigt viscoelastic model (parallel)

The same elements arranged in parallel give the
[Kelvin–Voigt model](https://en.wikipedia.org/wiki/Kelvin%E2%80%93Voigt_material):

```
     ---[viscous]---
 ----|              |----
     ---[elastic]---
```

```julia
kelvin_voigt = ParallelModel(viscous, elastic)
```

In a parallel model the *strain rate* is common to all elements, and the total *stress*
is the sum of contributions from each element.

## Example: visco-elasto-plastic model (series, three elements)

Adding a Drucker–Prager plastic element to the Maxwell model gives a VEP material:

```
---[viscous]---[elastic]---[plastic]---
```

```julia
viscous = LinearViscosity(1e22)
elastic = IncompressibleElasticity(10e9)
plastic = DruckerPrager(10e6, 30, 0)     # cohesion 10 MPa, ϕ=30°, ψ=0°

vep = SeriesModel(viscous, elastic, plastic)
```

See [`examples/Maxwell_VEP.jl`](https://github.com/albert-de-montserrat/RheologyCalculator.jl/blob/main/examples/Maxwell_VEP.jl)
for the full time-stepping example.

## Example: dilatant VEP model

When the dilatancy angle `ψ ≠ 0`, the plastic element also produces a volumetric strain
rate, which requires the fully compressible `Elasticity` type instead of
`IncompressibleElasticity`:

```julia
viscous = LinearViscosity(1e22)
elastic = Elasticity(10e9, 20e9)          # G = 10 GPa, K = 20 GPa
plastic = DruckerPrager(10e6, 30, 10)    # ψ = 10°

vep_dilatant = SeriesModel(viscous, elastic, plastic)
```

See [`examples/Maxwell_VEP_dilatant.jl`](https://github.com/albert-de-montserrat/RheologyCalculator.jl/blob/main/examples/Maxwell_VEP_dilatant.jl).

## Example: visco-elasto-viscoplastic model (VEVP, hybrid)

Elements can be nested to arbitrary depth.  The following represents a Maxwell backbone
with a regularised plastic element (plastic yield surface in parallel with a viscous
dashpot):

```
---[viscous]---[elastic]---+---[plastic]---+---
                           |               |
                           ---[viscous_reg]---
```

```julia
viscous     = LinearViscosity(1e22)
viscous_reg = LinearViscosity(1e20)     # regularisation viscosity
elastic     = IncompressibleElasticity(10e9)
plastic     = DruckerPrager(10e6, 0, 0)

p    = ParallelModel(plastic, viscous_reg)
vevp = SeriesModel(viscous, elastic, p)
```

See [`examples/Maxwell_VEVP.jl`](https://github.com/albert-de-montserrat/RheologyCalculator.jl/blob/main/examples/Maxwell_VEVP.jl).

## Example: Burgers viscoelastic model (hybrid)

The [Burgers material](https://en.wikipedia.org/wiki/Burgers_material) is a Maxwell
element in series with a Kelvin–Voigt element:

```
---[viscous₁]---[elastic₁]---+---[viscous₂]---+---
                              |                 |
                              ---[elastic₂]-----
```

```julia
viscous1 = LinearViscosity(1e21)
viscous2 = LinearViscosity(1e20)
elastic1 = Elasticity(1e10, 4e10)
elastic2 = Elasticity(1e10, 3e10)

p       = ParallelModel(viscous2, elastic2)
burgers = SeriesModel(viscous1, elastic1, p)
```

See [`examples/Burgers.jl`](https://github.com/albert-de-montserrat/RheologyCalculator.jl/blob/main/examples/Burgers.jl).

## The solver workflow

Once a composite model is assembled, the typical call sequence is:

```julia
# 1. Define constant (kinematic) variables
vars = (; ε = 1.0e-14, θ = 1.0e-20)   # deviatoric & volumetric strain rates

# 2. Define initial guesses for the differentiable unknowns
args = (; τ = 1e3, P = 1e6)

# 3. Non-differentiable auxiliary variables (updated each time step)
others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))

# 4. Build initial solution vector (warm start)
x = initial_guess_x(c, vars, args, others)

# 5. (Optional) build physics-based normalization vector
xnorm = normalisation_x(c, char_τ, char_ε)

# 6. Solve
x = solve(c, x, vars, others)           # without normalization
x = solve(c, x, vars, others, xnorm=xnorm)  # with normalization

# 7. Post-process: extract elastic stress/pressure for the next time step
τ_e = compute_stress_elastic(c, x, others)
P_e = compute_pressure_elastic(c, x, others)
```

### Variable roles

| Variable | Type | Role |
|:---|:---|:---|
| `vars` | `NamedTuple` | Kinematic inputs held **constant** during the local solve. Fields: `ε` (deviatoric strain rate, scalar or tuple), `θ` (volumetric strain rate). |
| `args` | `NamedTuple` | Initial guesses for the **differentiable unknowns** the solver iterates over (e.g. `τ`, `P`, `λ`). |
| `others` | `NamedTuple` | **Non-differentiable** auxiliary data needed by state functions: `dt` (time step), `τ0` (previous elastic deviatoric stress, tuple), `P0` (previous elastic pressure, tuple), `d` (grain size), etc. |
| `x` | `SVector` | Solution vector; one entry per equation in the composite model. |
| `xnorm` | `SVector` | Characteristic scales for each unknown, used to improve Newton convergence. |

### Normalization

For problems spanning many orders of magnitude (e.g. combining MPa stresses with
$\sim10^{-14}\ \mathrm{s}^{-1}$ strain rates), providing a normalization vector
significantly improves convergence:

```julia
char_τ  = plastic.C          # characteristic stress (e.g. cohesion)
char_ε  = vars.ε             # characteristic strain rate
xnorm   = normalisation_x(c, char_τ, char_ε)
x       = solve(c, x, vars, others, xnorm=xnorm)
```

`normalisation_x` assigns `char_τ` to stress/pressure unknowns and `char_ε` to
strain-rate/lambda unknowns based on the type of each equation.

