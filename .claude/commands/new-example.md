Create a new example script in `examples/` for RheologyCalculator.jl.

Ask the user:
1. What the example demonstrates (e.g. "Maxwell VEP with power-law viscosity", "Drucker-Prager cap under cyclic loading")
2. Which rheology elements to compose and how (series / parallel)
3. Whether to include a time-stepping loop or a single-point solve
4. Whether to plot results (if yes, use Makie — check existing examples for the plotting style)

Then create `examples/<DescriptiveName>.jl` following this structure:

```julia
using RheologyCalculator, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")   # lives in examples/

# -----------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------
# ... (material parameters, loading, timestep)

# -----------------------------------------------------------------------
# Composite model
# -----------------------------------------------------------------------
c = SeriesModel(...)   # or ParallelModel, or nested

# -----------------------------------------------------------------------
# Initial conditions
# -----------------------------------------------------------------------
vars   = (; ε = ...)
others = (; dt = ..., τ0 = 0.0, P0 = 0.0)
x      = initial_guess_x(c, vars, others)

# -----------------------------------------------------------------------
# Time loop
# -----------------------------------------------------------------------
nsteps = 500
τ_hist = zeros(nsteps)
t_hist = zeros(nsteps)

for i in 1:nsteps
    x = solve(c, x, vars, others; atol=1e-12)

    τ = x[1]
    τ_hist[i] = τ
    t_hist[i]  = i * others.dt

    # update elastic history
    τ0_new = compute_stress_elastic(c, x, vars, others)
    others  = merge(others, (; τ0 = τ0_new))
end

# -----------------------------------------------------------------------
# Plot (optional)
# -----------------------------------------------------------------------
# using GLMakie
# fig = Figure(); ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Stress (Pa)")
# lines!(ax, t_hist, τ_hist)
# display(fig)
```

Rules:
- Use `include("../rheologies/RheologyDefinitions.jl")` — not `using` — so the script is self-contained
- Always update `τ0` (and `P0` if volumetric) via `compute_stress_elastic` / `compute_pressure_elastic` each timestep
- Use `atol=1e-12` in `solve` unless there's a specific reason to relax it
- Do not add `Pkg.activate` or environment setup — the user handles that
- Keep parameters at the top in a clearly labeled block, not scattered through the script
- No file I/O or hardcoded paths
- If plotting: use GLMakie, follow the style in `examples/Maxwell_VEP.jl`
- Name the file descriptively in PascalCase matching the pattern of existing examples
