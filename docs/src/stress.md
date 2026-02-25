# Usage examples

The following examples are drawn from the `examples/` folder and illustrate the
complete time-stepping workflow for different composite rheologies.

---

## Maxwell viscoelastic model

**File:** `examples/Maxwell_VE.jl`

The [Maxwell model](https://en.wikipedia.org/wiki/Maxwell_material) couples a viscous
damper and an elastic spring in **series**.  The analytical stress relaxation for a
constant strain rate $\dot\varepsilon$ is:

$$\tau(t) = 2\dot\varepsilon\,\eta \left(1 - e^{-Gt/\eta}\right)$$

### Setup

```julia
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic
include("../rheologies/RheologyDefinitions.jl")

viscous = LinearViscosity(1e22)
elastic = IncompressibleElasticity(10e9)
c       = SeriesModel(viscous, elastic)

vars   = (; ε = 1.0e-14, θ = 1.0e-20)              # constant deviatoric & volumetric strain rate
args   = (; τ = 2.0e3, P = 1.0e6)                  # initial guess for τ and P
others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,)) # Δt and previous elastic state

x = initial_guess_x(c, vars, args, others)
```

### Time loop

```julia
function stress_time(c, vars, x; ntime = 1000, dt = 1.0e10)
    τ_history = zeros(ntime)
    t_history = zeros(ntime)
    τ_e = (0.0,)  # previous elastic deviatoric stress (one entry per elastic element)
    t   = 0.0

    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, P0 = (0.0,))

        x      = solve(c, x, vars, others)      # Newton-Raphson solve
        τ_e    = compute_stress_elastic(c, x, others)  # update elastic state

        τ_history[i] = x[1]
        t            += dt
        t_history[i]  = t
    end
    return t_history, τ_history
end

t_v, τ = stress_time(c, vars, x)
```

`τ0` is a **tuple** with one entry for each elastic element in the model; here there is
one `IncompressibleElasticity`, so `τ0 = (τ_xx_elastic,)`.

---

## Maxwell visco-elasto-plastic model (VEP)

**File:** `examples/Maxwell_VEP.jl`

Adds a Drucker-Prager plastic element in **series** to the Maxwell model.  The stress
is capped at the yield surface once the material reaches plastic failure.

```julia
viscous = LinearViscosity(1e22)
elastic = IncompressibleElasticity(10e9)
plastic = DruckerPrager(10e6, 30, 0)   # cohesion 10 MPa, ϕ=30°, ψ=0° (non-dilatant)

c = SeriesModel(viscous, elastic, plastic)

vars   = (; ε = 1.0e-14, θ = 1.0e-20)
args   = (; τ = 0.0, λ = 0.0)                       # λ is the plastic multiplier
others = (; dt = 1.0e8, P = 1.0e6, τ0 = 0.0, P0 = 0.0)

x      = initial_guess_x(c, vars, args, others)

# Use a normalization vector for better Newton convergence
char_τ = plastic.C
char_ε = vars.ε
xnorm  = normalisation_x(c, char_τ, char_ε)

x = solve(c, x, vars, others, xnorm=xnorm)
```

---

## Dilatant VEP model

**File:** `examples/Maxwell_VEP_dilatant.jl`

When `ψ ≠ 0`, plastic flow generates a volumetric strain rate. This requires:

- `Elasticity(G, K)` instead of `IncompressibleElasticity(G)` (to handle volumetric elasticity).
- `θ` (volumetric strain rate) in `vars`.
- `P0` updated each step via `compute_pressure_elastic`.

```julia
viscous = LinearViscosity(1e22)
elastic = Elasticity(10e9, 20e9)
plastic = DruckerPrager(10e6, 30, 10)   # ψ = 10° → dilatant

c = SeriesModel(viscous, elastic, plastic)

vars   = (; ε = 1.0e-14, θ = 1.0e-20)
args   = (; τ = 0.0, P = 1.0e6, λ = 0.0)
others = (; dt = 1.0e8, τ0 = 0.0, P0 = 0.0)

x = initial_guess_x(c, vars, args, others)
```

Inside the time loop, both the elastic stress **and pressure** must be updated:
```julia
τ_e = compute_stress_elastic(c, x, others)
P_e = compute_pressure_elastic(c, x, others)
others = (; dt = dt, τ0 = τ_e, P0 = P_e)
```

---

## Visco-elasto-viscoplastic model (VEVP)

**File:** `examples/Maxwell_VEVP.jl`

A **regularised** plastic element — plasticity in *parallel* with a small viscosity —
prevents the singular yield surface from causing solver instability:

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
c    = SeriesModel(viscous, elastic, p)

vars   = (; ε = 1.0e-14)
args   = (; τ = 2.0e3, λ = 0.0)
others = (; dt = 1.0e8, τ0 = (0.0,), P = 1.0e6, P0 = (0.0,))

x     = initial_guess_x(c, vars, args, others)
xnorm = normalisation_x(c, 1e6, vars.ε)

x = solve(c, x, vars, others, xnorm=xnorm)
```

---

## Burgers viscoelastic model

**File:** `examples/Burgers.jl`

The [Burgers model](https://en.wikipedia.org/wiki/Burgers_material) is a Maxwell element
in series with a Kelvin-Voigt element:

```
---[viscous₁]---[elastic₁]---+---[viscous₂]---+---
                              |                 |
                              ---[elastic₂]-----
```

This model has **two** elastic elements, so `τ0` and `P0` in `others` are tuples of
length 2:

```julia
viscous1 = LinearViscosity(1e21)
viscous2 = LinearViscosity(1e20)
elastic1 = Elasticity(1e10, 4e10)
elastic2 = Elasticity(1e10, 3e10)

p       = ParallelModel(viscous2, elastic2)
burgers = SeriesModel(viscous1, elastic1, p)

vars   = (; ε = 1.0e-15, θ = 1.0e-20)
args   = (; τ = 2.0e3, P = 1.0e6)
others = (; dt = 1.0e10, τ0 = (1.0, 2.0), P0 = (0.0, 0.1))

x = initial_guess_x(burgers, vars, args, others)
```

Inside the time loop:
```julia
τ_e = compute_stress_elastic(burgers, x, others)   # returns a tuple of length 2
P_e = compute_pressure_elastic(burgers, x, others) # returns a tuple of length 2
others = (; dt = dt, τ0 = τ_e, P0 = P_e)
```

---

## VEP + Cap model (Popov et al., 2025)

**File:** `examples/Maxwell_VEPCap.jl`

Implements the mode-1/mode-2 plasticity model with a tensile cap (Popov et al., 2025,
*EGUsphere*).  The `DruckerPragerCap` rheology is defined in
`rheologies/DruckerPragerCap.jl`:

```julia
include("../rheologies/DruckerPragerCap.jl")

viscous = LinearViscosity(1e23)
elastic = Elasticity(1e10, 2e11)
plastic = DruckerPragerCap(; C=1e6, ϕ=30.0, ψ=10.0, η_vp=0.0, Pt=-5e5)

c = SeriesModel(viscous, elastic, plastic)

vars   = (; ε = 0.0, θ = 7.0e-15)
args   = (; τ = 0.0, P = 0.3e6, λ = 0.0)
others = (; dt = 1.0e5, τ0 = (0.0,), P0 = (0.3e6,))

x     = initial_guess_x(c, vars, args, others)
xnorm = normalisation_x(c, plastic.C, vars.ε + vars.θ)
```

---

## Solver API summary

### `initial_guess_x`

```julia
x = initial_guess_x(c, vars, args, others)
```

Constructs a physics-based warm-start `SVector` for the Newton solver.  Stress unknowns
are estimated as an arithmetic sum of individual element responses; strain-rate unknowns
as a harmonic mean.

### `normalisation_x`

```julia
xnorm = normalisation_x(c, char_τ, char_ε)
```

Returns a `SVector` of characteristic scales for each unknown:
`char_τ` for stress/pressure entries, `char_ε` for strain-rate/lambda entries.
Providing `xnorm` to `solve` via the `xnorm` keyword argument improves Newton
convergence for problems spanning many orders of magnitude.

### `solve`

```julia
x = solve(c, x, vars, others;
          xnorm   = nothing,   # optional normalization vector
          atol    = 1.0e-9,    # absolute residual tolerance
          rtol    = 1.0e-9,    # relative residual tolerance
          itermax = 1.0e4,     # maximum Newton iterations
          verbose = false)     # print iteration info
```

Newton-Raphson solver with automatic Jacobian (ForwardDiff) and backtracking line
search.  Returns the converged solution vector `x`.

### `compute_stress_elastic` / `compute_pressure_elastic`

```julia
τ_e = compute_stress_elastic(c, x, others)
P_e = compute_pressure_elastic(c, x, others)
```

Post-processing functions that extract the elastic deviatoric stress and pressure from
the converged solution `x`.  The returned tuples have one entry per elastic element in
the composite model; they should be stored and passed as `τ0` / `P0` in `others` at
the next time step.
