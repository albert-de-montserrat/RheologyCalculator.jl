Add a test file for a rheology element or composite model in RheologyCalculator.jl.

Ask the user:
1. Which rheology or composite to test (e.g. `DruckerPrager`, `Maxwell VEP`)
2. Whether an analytical solution exists to compare against
3. Whether to test series, parallel, or both compositions
4. Typical parameter values

Then create `test/test_<Name>.jl` following this structure:

```julia
using RheologyCalculator, Test, Statistics
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import RheologyCalculator: compute_residual

include("../rheologies/RheologyDefinitions.jl")
include("../examples/tensor_helpers.jl")

@testset "<Name>" begin

    # --- Parameters ---
    # ... (fill in)

    # --- Build composite ---
    c = SeriesModel(...)

    # --- Initial conditions ---
    vars   = (; ε = ...)
    others = (; dt = ..., τ0 = 0.0, P0 = 0.0)
    x      = initial_guess_x(c, vars, others)

    # --- Time stepping ---
    nsteps = 100
    τ_hist = zeros(nsteps)
    for i in 1:nsteps
        x = solve(c, x, vars, others)
        τ = x[1]
        τ_hist[i] = τ
        # update history
        others = merge(others, (; τ0 = τ))
    end

    # --- Assertions ---
    # Residual should be near zero at final state
    r = compute_residual(c, x, vars, others)
    @test maximum(abs, r) < 1e-10

    # If analytical solution known, compare:
    # τ_analytical = ...
    # @test isapprox(τ_hist[end], τ_analytical, rtol=1e-6)

end
```

Rules:
- The file will be auto-discovered by `test/runtests.jl` — no registration needed
- Always test residual norm at the final state
- For plastic models, include a test that λ ≥ 0 (plastic multiplier non-negative)
- For elastic models, include a creep-recovery or stress-relaxation check
- For volumetric models, also check pressure residuals
- Keep `nsteps` small (≤ 200) so the test suite stays fast
- Do not use `@test_broken` — if a case isn't implemented, omit it

After creating the file, remind the user to run:
```
julia --project=. -e 'include("test/runtests.jl")'
```
or just `include("test/test_<Name>.jl")` in the REPL to verify it passes.
