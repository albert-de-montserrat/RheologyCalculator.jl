# SKILLS — RheologyCalculator.jl

Project-specific patterns and workflows for working in this codebase.

---

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Individual test files can be included directly in the REPL after activating the project:

```julia
using Pkg; Pkg.activate(".")
include("test/runtests.jl")   # runs all test_*.jl files
include("test/test_VEP.jl")   # run a single suite
```

The test runner (`test/runtests.jl`) auto-discovers every file prefixed `test_` in the test directory.

---

## Package structure

| Path | Purpose |
|---|---|
| `src/rheology_types.jl` | Abstract type hierarchy (`AbstractRheology`, `AbstractViscosity`, `AbstractElasticity`, `AbstractPlasticity`) and state-function dispatch machinery |
| `src/composite.jl` | `SeriesModel` / `ParallelModel` — leaf/branch separation and state-function routing |
| `src/equations.jl` | `generate_equations` + `compute_residual` — builds and evaluates the nonlinear residual system |
| `src/solver.jl` | Newton-Raphson solver with backtracking line search (`solve`) |
| `src/state_functions.jl` | Canonical state-function names and their fallback/wrapper boilerplate |
| `src/initial_guess.jl` | `initial_guess_x`, `x_keys` |
| `src/normalize_x.jl` | `normalisation_x` |
| `src/strain_rate_correction.jl` | `effective_strain_rate_correction` |
| `rheologies/RheologyDefinitions.jl` | Concrete rheology elements (template for user-defined laws) |
| `rheologies/` | Additional plasticity models (Drucker-Prager Cap, ModCamClay, Hyperbolic, Golchin, RateState) |
| `examples/` | End-to-end driver scripts (Maxwell VE, VEP, VEVP, dislocation creep, etc.) |
| `test/` | Test suites per rheology type |

---

## Defining a new rheology element

1. **Subtype** one of `AbstractViscosity`, `AbstractElasticity`, or `AbstractPlasticity`.
2. **Declare state functions** — return which canonical functions this element contributes to:
   ```julia
   @inline series_state_functions(::MyRheology)   = (compute_strain_rate,)
   @inline parallel_state_functions(::MyRheology) = (compute_stress,)
   ```
3. **Implement** the corresponding function(s):
   ```julia
   @inline compute_strain_rate(r::MyRheology; τ = 0, kwargs...) = τ / (2r.η)
   @inline compute_stress(r::MyRheology; ε = 0, kwargs...) = 2r.η * ε
   ```
4. For **elastic** rheologies, also implement `compute_viscosity` and handle history kwargs (`τ0`, `dt`).
5. For **plastic** rheologies, implement `compute_lambda` (series) or `compute_lambda_parallel` (parallel) for the yield consistency equation, plus `compute_plastic_strain_rate`.
6. Declare `history_kwargs` if the element reads element-local history fields from `others` (e.g., `τ0` per element):
   ```julia
   @inline history_kwargs(::MyRheology) = (:τ0,)
   ```
7. For **volumetric** behavior, add `_isvolumetric(::MyRheology) = Val(true)` and implement `compute_volumetric_strain_rate` / `compute_pressure`.

See `rheologies/RheologyDefinitions.jl` for concrete examples covering all element types.

---

## Building and solving a composite model

```julia
include("rheologies/RheologyDefinitions.jl")

# Series Maxwell VE
c = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))

# Initial solver vector
x = initial_guess_x(c, vars, others)

# Solve
vars  = (; ε = 1e-15)                          # prescribed (differentiable) inputs
others = (; dt = 1e10, τ0 = 0.0, P0 = 0.0)   # non-differentiable auxiliary state

x = solve(c, x, vars, others; atol=1e-12, verbose=true)
```

- `vars` holds prescribed inputs that become the right-hand side of the global equations (e.g., total strain rate `ε`, volumetric strain rate `θ`).
- `others` holds non-differentiable fields (timestep `dt`, previous-step stress `τ0`, pressure `P0`, material parameters used only as numbers, not as AD seeds).
- The solver uses **ForwardDiff** to compute the Jacobian; all residual code must be AD-compatible.
- For tensor inputs, pass `ε` as an `NTuple{3}` or `NTuple{6}` (Voigt notation); the solver reduces to the second invariant internally via `effective_strain_rate_correction`.

---

## Adding a test

Create `test/test_MyRheology.jl`. The runner picks it up automatically. Typical structure:

```julia
using RheologyCalculator, Test

include("../rheologies/RheologyDefinitions.jl")
include("../examples/tensor_helpers.jl")

@testset "MyRheology" begin
    c   = SeriesModel(...)
    x   = initial_guess_x(c, vars, others)
    x   = solve(c, x, vars, others)
    # assert stress / strain-rate consistency
    @test ...
end
```

---

## Key design constraints

- **No heap allocations in the hot path.** The equation-generation loop and residual evaluation are fully `@generated` / `@inline` and work on `SVector` / `NTuple` types. Do not introduce heap-allocating containers in solver-critical code.
- **ForwardDiff compatibility.** Every function called inside `compute_residual` must propagate dual numbers. Avoid `if/else` on numeric values; prefer `ifelse` or smooth formulations.
- **State-function dispatch is compile-time.** The set of equations is determined purely from the type of the composite model — there is no runtime branching over element types.
- **`vars` vs `others` separation matters.** Only fields in `vars` appear in the Jacobian columns. Fields that should not be differentiated (geometry, history from the previous timestep) must go into `others`.
- **Element-local history** is dispatched via `history_kwargs` + `extract_local_kwargs`. When a composite contains multiple elastic or plastic elements, pass their history fields as tuples: `others = (; τ0 = (τ0_el1, τ0_el2), dt = 1e10)`.

---

## Plasticity models available

| File | Model |
|---|---|
| `rheologies/RheologyDefinitions.jl` | `DruckerPrager`, `VonMises` |
| `rheologies/DruckerPragerCap.jl` | `DruckerPragerCap` (Popov et al., 2025) |
| `rheologies/Hyperbolic.jl` | Hyperbolic smoothed Drucker-Prager (Abbo & Sloan, 1995) |
| `rheologies/ModCamClay.jl` | Modified Cam-Clay (classical + Golchin variant) |
| `rheologies/RateState_HypoPlastic.jl` | Rate-and-State friction (Herrendörfer et al., 2018) |

---

## Debugging tips

- Use `solve_timed(c, x, vars, others)` (`src/solver.jl`) to get a `TimerOutputs` breakdown of Newton-Raphson overhead.
- Use `generate_equations(c)` to inspect the equation tree without running the solver.
- Use `compute_residual(c, x, vars, others)` directly to check residuals at any `x`.
- Print the composite with `show(c)` — `print_rheology.jl` provides a structured display.
