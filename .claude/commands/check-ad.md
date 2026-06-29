Verify that a rheology element's state functions are ForwardDiff-compatible.

Ask the user which rheology to check (e.g. `DruckerPrager`, `LinearViscosity`).

Then:

1. Read the rheology's implementation in `rheologies/RheologyDefinitions.jl` (or the relevant file in `rheologies/`).

2. Identify every state function it implements (`compute_strain_rate`, `compute_stress`, `compute_lambda`, `compute_viscosity`, etc.).

3. For each function, check for these common AD-breaking patterns:
   - `if condition` or `cond ? a : b` where `condition` involves a numeric value that will be a dual number (use `ifelse` instead)
   - `abs(x)` in a context where `x` can be zero (subgradient issue — use `sqrt(x^2 + ε)` or `hypot`)
   - `sign(x)` — not differentiable at zero; flag if used without a smooth workaround
   - `convert`, `Int(...)`, `round`, `floor`, `ceil` applied to a differentiable variable
   - `@assert` or `error(...)` inside a hot path (breaks tracing)
   - Mutation of arrays (`.=`, `push!`, etc.) — incompatible with ForwardDiff dual arrays
   - Calls to external functions that are not known to be AD-compatible

4. Write a small Julia snippet the user can run in the REPL to empirically verify AD compatibility:

```julia
using ForwardDiff
include("rheologies/RheologyDefinitions.jl")

r = MyRheology(...)   # fill in typical parameter values
τ_test = 1e4          # typical stress value

# Test compute_strain_rate
f = τ -> compute_strain_rate(r; τ = τ)
println("dε/dτ = ", ForwardDiff.derivative(f, τ_test))

# Add similar snippets for compute_stress, compute_lambda, etc.
```

5. Report findings:
   - List any problematic patterns found (file:line)
   - Confirm which functions passed the pattern check
   - Note if the empirical snippet is expected to pass or likely to error

Do not run the snippet yourself — output it for the user to paste into the REPL.
