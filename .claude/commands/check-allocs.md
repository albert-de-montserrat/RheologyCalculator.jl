Check for heap allocations in the hot path of a RheologyCalculator.jl composite model.

Ask the user:
1. Which composite model to benchmark (or offer to use an example from `examples/`)
2. Typical parameter values to use (viscosity, elastic modulus, timestep, strain rate, etc.)

Then write a self-contained Julia benchmark script the user can run in the REPL:

```julia
using RheologyCalculator, StaticArrays, BenchmarkTools
include("rheologies/RheologyDefinitions.jl")

# --- User fills in their composite ---
c      = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10, 1e10))
vars   = (; ε = 1e-15)
others = (; dt = 1e10, τ0 = 0.0, P0 = 0.0)
x      = initial_guess_x(c, vars, others)

# Warm up (force compilation)
x = solve(c, x, vars, others)

# Allocation check
allocs = @allocated solve(c, x, vars, others)
println("solve allocations: ", allocs, " bytes")

# Residual-only check
allocs_r = @allocated RheologyCalculator.compute_residual(c, x, vars, others)
println("compute_residual allocations: ", allocs_r, " bytes")

# Full benchmark
@btime solve($c, $x, $vars, $others)
```

After writing the script, inspect the composite model type the user described and flag any likely allocation sources:
- Non-`@generated` functions over variable-length tuples
- `NTuple` operations that fall back to generic iteration
- `NamedTuple` merges inside generated code
- Any `Vector` or `Dict` usage in state functions or kwargs helpers
- `compute_residual` returning `SA[residual...]` — check whether `residual` is a fixed-length `NTuple` (good) or dynamic (bad)

Report which functions are suspects and suggest fixes if allocations are found (typically: ensure the composite type is fully concrete so `generate_equations` can specialize, or add `@inline` / `@generated` to the offending method).

Do not run the benchmark yourself — output it for the user to paste into the REPL.
