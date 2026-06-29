Identify uncovered code in RheologyCalculator.jl and write targeted tests to improve coverage.

## Step 1 — Run the test suite with coverage

```bash
julia --project=. --code-coverage=user --startup-file=no test/runtests.jl
```

This writes `.cov` files next to each source file in `src/`. Do NOT use
`Pkg.test()` — it precompiles with `--code-coverage=none` and generates no
`.cov` files. If the user says they already have fresh `.cov` files, skip this
step.

## Step 2 — Parse coverage data

Run the following snippet in the REPL (or use the Bash tool) to find uncovered lines:

```julia
using Coverage
coverage = process_folder("src")
covered, total = get_summary(coverage)
println("Coverage: $(round(100*covered/total, digits=1))% ($covered / $total lines)")

# Print uncovered lines grouped by file
for fc in coverage
    uncov = [l.line for l in fc.coverage if l !== nothing && l.hits == 0]
    isempty(uncov) && continue
    println("\n$(fc.filename):")
    for ln in uncov
        println("  line $ln")
    end
end
```

If `Coverage.jl` is not in the project, add it temporarily:
```bash
julia --project=. -e 'using Pkg; Pkg.add("Coverage")'
```

## Step 3 — Identify meaningful gaps

Read each uncovered file at the flagged lines. Skip:
- Lines inside `@generated` function bodies (these are macro-expansion internals, not runtime paths)
- `@inline` fallbacks that intentionally dead-end (`error(...)`, `throw(...)`)
- Lines in files that are not `include`d by `src/RheologyCalculator.jl`

Focus on:
- Functions reachable from user-facing APIs (`solve`, `initial_guess_x`, `normalisation_x`, `compute_stress_elastic`, `compute_pressure_elastic`, `effective_strain_rate_correction`)
- Branches that differ by composition topology (e.g. `ParallelModel` containing a `SeriesModel` branch, or a model with `isvolumetric = true`)
- `kwargs.jl` helpers (`split_args`, `augment_args`, `all_differentiable_kwargs`)
- Normalisation dispatch for plastic unknowns (`:λ`, `:τ_pl`)
- `post_calculations.jl` elastic history extraction for stress-form vs. strain-rate-form equations

## Step 4 — Write new tests

Create or extend test files under `test/`. Each new test file must be named `test_<Topic>.jl`
so `test/runtests.jl` picks it up automatically.

Standard file header:
```julia
using RheologyCalculator, Test, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import RheologyCalculator: compute_residual

include("../rheologies/RheologyDefinitions.jl")
include("../examples/tensor_helpers.jl")
```

### Coverage checklist per gap type

**`kwargs.jl` — `split_args`, `differentiable_kwargs`, `normalisation_x` for plastic unknowns**
```julia
@testset "kwargs helpers" begin
    c    = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10), DruckerPrager(10e6, 30, 0))
    eqs  = RheologyCalculator.generate_equations(c)
    # normalisation covers :λ (plastic multiplier) → char_τ path
    xn   = normalisation_x(c, 1e6, 1e-15)
    @test length(xn) == length(eqs)
    # x_keys covers all canonical unknowns
    @test :λ in x_keys(c)
end
```

**`others.jl` — `isvolumetric` for composites with volumetric elements**
```julia
@testset "isvolumetric" begin
    c_inc = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
    c_vol = SeriesModel(LinearViscosity(1e22), Elasticity(1e10, 1e10))
    @test RheologyCalculator.isvolumetric(c_inc) == Val(false)
    @test RheologyCalculator.isvolumetric(c_vol) == Val(true)
end
```

**`post_calculations.jl` — `compute_stress_elastic` and `compute_pressure_elastic`**
```julia
@testset "elastic history extraction" begin
    c      = SeriesModel(LinearViscosity(1e22), Elasticity(1e10, 1e10))
    vars   = (; ε = 1e-15, θ = 1e-20)
    args   = (; τ = 0.0, P = 0.0)
    others = (; dt = 1e10, τ0 = (0.0,), P0 = (0.0,))
    x      = initial_guess_x(c, vars, args, others)
    x      = solve(c, x, vars, others)

    τ_el = compute_stress_elastic(c, x, others)
    P_el = compute_pressure_elastic(c, x, others)
    @test length(τ_el) == 1
    @test length(P_el) == 1
    @test all(isfinite, τ_el)
    @test all(isfinite, P_el)
end
```

**`strain_rate_correction.jl` — tensor NTuple path of `effective_strain_rate_correction`**
```julia
@testset "strain rate correction (tensor)" begin
    c      = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
    ε      = (1e-15, -1e-15, 0.0)          # 2-D Voigt tensor
    others = (; dt = 1e10, τ0 = (0.0,), P0 = (0.0,))
    δε     = effective_strain_rate_correction(c, ε, others.τ0, others)
    @test length(δε) == length(ε)
    @test all(isfinite, δε)
end
```

**`initial_guess.jl` — `tensor2invariant` for nested NTuple and NamedTuple inputs**
```julia
@testset "tensor2invariant" begin
    import RheologyCalculator: tensor2invariant
    @test tensor2invariant(1.0) ≈ 1.0
    @test tensor2invariant((1e-15, -1e-15, 0.0)) isa Float64
    nt = (; ε = (1e-15, -1e-15, 0.0), τ = 1e3)
    nt2 = tensor2invariant(nt)
    @test nt2.τ ≈ 1e3          # scalar passes through
    @test nt2.ε isa Float64    # tensor reduced to invariant
end
```

## Step 5 — Remove .cov files

Delete all `.cov` files generated by the coverage run. The test runner
(`test/runtests.jl`) auto-discovers files prefixed `test_`, and `.cov` files
in `test/` match that prefix, causing load errors if left behind.

```bash
find . -name "*.cov" -delete
```

## Step 6 — Verify

Run the new tests and confirm the suite passes:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Then re-run with `--code-coverage=user` and check that the coverage percentage increased
relative to the baseline from Step 2.

## Rules

- Only write tests for code that is actually `include`d in `src/RheologyCalculator.jl`
- Keep each test small and self-contained; use `let` blocks to avoid polluting global scope
- Always assert `isfinite` on solver outputs — silent `NaN` propagation is the most common silent failure
- Do not test private `_`-prefixed helpers directly; reach them through a public API call
- Do not aim for 100% line coverage of `@generated` expansions — that's not feasible without deep metaprogramming introspection
