# RheologyCalculator.jl

`RheologyCalculator.jl` is a Julia package for computing the **stress-strain response of
arbitrary composite rheological materials** at a single integration point. It is designed
for use inside larger finite-element or finite-difference codes (e.g. geodynamic solvers)
where the local constitutive update must be solved cheaply at every quadrature point.

## Features

- **Composable rheology** — assemble any combination of viscous, elastic, and plastic
  elements in series, in parallel, or in arbitrarily nested hybrid configurations.
- **Automatic Newton–Raphson solver** — the system of nonlinear equations is solved with a
  backtracking Newton–Raphson method with automatic Jacobian evaluation via
  [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).
- **Allocation-free inner loop** — solution vectors are stored as `SVector`s from
  [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl); internal loops are
  fully unrolled at compile time through `@generated` functions, making the kernel
  suitable for GPU execution or tight CPU loops.
- **Extensible rheology** — adding a new creep law or plasticity model requires only
  defining a struct and the relevant `compute_*` state functions; the equation-generation
  and solver machinery is inherited automatically.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/albert-de-montserrat/RheologyCalculator.jl")
```

## Quick-start: Maxwell viscoelastic model

Below is the minimal workflow to compute the stress evolution of a **Maxwell
viscoelastic** model (viscous element in series with an incompressible elastic spring).

```julia
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic

include("rheologies/RheologyDefinitions.jl")  # defines individual rheology types

# 1. Define individual rheological elements
viscous = LinearViscosity(1e22)            # η = 1e22 Pa·s
elastic = IncompressibleElasticity(10e9)   # G = 10 GPa

# 2. Build the composite model (series = Maxwell)
c = SeriesModel(viscous, elastic)

# 3. Define input (constant) and unknown (differentiable) variables
vars   = (; ε = 1.0e-14, θ = 1.0e-20)              # deviatoric & volumetric strain rates
args   = (; τ = 2.0e3, P = 1.0e6)                  # initial guesses for τ and P
others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)) # time step and previous elastic stresses

# 4. Compute initial guess for the solution vector
x = initial_guess_x(c, vars, args, others)

# 5. Solve for τ
x = solve(c, x, vars, others)

# 6. Extract elastic stress (needed to carry forward in time)
τ_e = compute_stress_elastic(c, x, others)
```

For a full time-stepping loop and a comparison against the analytical Maxwell solution,
see [Stress–time curve of a Burger's material](stress.md) or the
[`examples/Maxwell_VE.jl`](https://github.com/albert-de-montserrat/RheologyCalculator.jl/blob/main/examples/Maxwell_VE.jl)
script.

## Package layout

| Module path | Purpose |
|:---|:---|
| `src/composite.jl` | `SeriesModel`, `ParallelModel` composite types |
| `src/rheology_types.jl` | Abstract type hierarchy (`AbstractRheology`, subtypes) |
| `src/state_functions.jl` | `compute_strain_rate`, `compute_stress`, … fallback dispatch |
| `src/equations.jl` | Automatic equation generation from a composite model |
| `src/initial_guess.jl` | `initial_guess_x` — warm-starting the Newton solver |
| `src/normalize_x.jl` | `normalisation_x` — physics-based scaling of the solution vector |
| `src/solver.jl` | Newton–Raphson solver with backtracking line search |
| `src/post_calculations.jl` | `compute_stress_elastic`, `compute_pressure_elastic` |
| `rheologies/RheologyDefinitions.jl` | Reference implementations of all built-in creep laws |

## Contents

```@contents
Pages = ["composites.md", "rheology.md", "stress.md"]
Depth = 2
```
