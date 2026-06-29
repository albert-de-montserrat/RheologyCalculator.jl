Audit and write Julia docstrings for RheologyCalculator.jl functions, types, and methods.

Ask the user:
1. Which file(s) or symbol(s) to document (e.g. "all of src/solver.jl", "DruckerPrager", "compute_residual")
2. Whether to audit existing docstrings for completeness or write new ones from scratch

## Docstring style for this project

Follow the Julia docstring conventions used in this codebase (see `src/solver.jl`, `src/equations.jl`, `src/state_functions.jl` for examples):

```julia
"""
    function_name(arg1, arg2; kwarg=default)

One-line summary of what the function returns or does.

Longer explanation only when the behaviour is non-obvious. Explain the
*why* and any constraints, not what the code does line by line.

# Arguments / Fields
- `arg1`: description
- `kwarg`: description (include default and when to change it)

# Example
```julia
c = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
x = solve(c, x, vars, others)
```
"""
```

Rules:
- One-liner is mandatory. Everything else is optional.
- Do NOT explain what the code does — explain what it *means* and any non-obvious invariants or constraints.
- For structs: document each field only when the name isn't self-explanatory.
- For state functions (`compute_strain_rate`, `compute_lambda`, etc.): explain which composition context uses this (series vs. parallel) and what the solver does with the return value.
- For `@generated` or `@inline` methods: document the *non-generated* conceptual version; implementation details of the generated expansion are not useful in docs.
- Keep examples minimal — one or two lines showing the typical call site.
- Use `@ref` cross-links when mentioning other documented symbols: `[`SeriesModel`](@ref)`.
- Do not add docstrings to private helpers (names starting with `_`) unless they are exported or referenced from `docs/src/api.md`.

## What to check in an audit

When auditing existing docstrings, flag:
- Missing docstring on an exported symbol
- Docstring that describes *what* the code does rather than *why* or *what it means*
- Signature in the docstring that doesn't match the current method signature
- Missing `# Example` on a user-facing function (solver, composite constructors, state functions)
- `@docs` block in `docs/src/api.md` referencing a symbol that has no docstring (will cause Documenter to warn)

## Documenter integration

After writing docstrings, check that the symbol is listed in the appropriate `@docs` block in `docs/src/`:
- Public API → `docs/src/api.md`
- Composite types → `docs/src/composites.md`
- Rheology types and interface → `docs/src/rheology.md`
- Solver and post-processing → `docs/src/stress.md`

If a new symbol should appear in the docs but isn't listed, add it to the correct `@docs` block. `make.jl` uses `checkdocs = :exports` so every exported symbol needs a docstring or Documenter will warn.
