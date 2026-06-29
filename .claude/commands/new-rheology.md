Scaffold a new rheology element for RheologyCalculator.jl.

Ask the user for:
1. Element name (e.g. `DruckerPrager`, `PowerLawViscosity`)
2. Type: viscosity / elasticity / plasticity
3. Parameters (names and types, e.g. `η::T`, `G::T, dt::T`)
4. Is it volumetric? (needs pressure/volumetric strain-rate equations)
5. Does it carry history fields? (e.g. `τ0`, `P0` from the previous timestep)
6. Series or parallel or both?

Then generate the following in `rheologies/RheologyDefinitions.jl` (or a new file in `rheologies/` if it's a standalone model):

- A `struct` subtyping `AbstractViscosity`, `AbstractElasticity`, or `AbstractPlasticity`
- `@inline series_state_functions(::MyRheology) = (...)` — include `compute_strain_rate` for viscous/elastic; for plastic add `compute_lambda`
- `@inline parallel_state_functions(::MyRheology) = (...)` — include `compute_stress` for viscous/elastic; for plastic add `compute_lambda_parallel` and `compute_plastic_strain_rate`
- All state-function methods using keyword-argument style: `compute_strain_rate(r::MyRheology; τ=0, kwargs...) = ...`
- `compute_viscosity`, `compute_viscosity_series`, `compute_viscosity_parallel` for viscous/elastic elements
- If elastic: `history_kwargs(::MyRheology) = (:τ0,)` (add `:P0` if volumetric); use `τ0` and `dt` inside `compute_strain_rate`
- If volumetric: `_isvolumetric(::MyRheology) = Val(true)` plus `compute_volumetric_strain_rate` and `compute_pressure`
- If plastic (series): `compute_lambda` returning the yield-function residual
- If plastic (parallel): `compute_lambda_parallel` and `compute_plastic_strain_rate`

Rules to follow:
- All state functions must be `@inline` and AD-compatible (no `if/else` on numeric values; use `ifelse` or smooth approximations)
- Keyword methods must accept `kwargs...` and ignore unknown keys
- Import any new function names at the top of the file with `import RheologyCalculator: ...`
- Do not add heap-allocating containers
- Match the code style in `rheologies/RheologyDefinitions.jl` (no comments explaining what the code does, only non-obvious WHY comments)

After writing the code, remind the user to add a corresponding test in `test/test_<Name>.jl` using `/add-test`.
