Build, update, or extend the Documenter.jl documentation for RheologyCalculator.jl.

## Documentation structure

```
docs/
  make.jl          ← Documenter entry point; lists pages and modules
  src/
    index.md       ← Landing page / overview
    composites.md  ← SeriesModel, ParallelModel, CompositeModel
    rheology.md    ← Abstract types, element interface, concrete creep/plastic laws
    stress.md      ← Solver, initial guess, normalisation, post-processing
    api.md         ← Full @docs autodoc blocks for all exported and internal symbols
  assets/          ← PNG figures referenced in pages
```

## How to build locally

```bash
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html` in a browser. The build uses `warnonly = Documenter.except(:footnote)` so most warnings are fatal — fix them before pushing.

## Common tasks

### Add a new page

1. Create `docs/src/mypage.md`.
2. Add `"My Page" => "mypage.md"` to the `pages` vector in `docs/make.jl`.
3. Populate the page with prose and `@docs` blocks as needed.

### Add a new symbol to the API page

Open `docs/src/api.md` and add the fully-qualified symbol to the appropriate `@docs` block:

```markdown
```@docs
RheologyCalculator.MyNewFunction
```
```

Use `RheologyCalculator.Symbol` for internal (non-exported) symbols; bare `Symbol` for exported ones.

### Add a new rheology element to the Rheology page

In `docs/src/rheology.md`, add a subsection under the appropriate category (Creep Laws / Elasticity / Plastic Failure) with:
- The constructor signature in a Julia code block
- The governing equation in a `math` block (LaTeX)
- A one-line description of parameters

### Fix a Documenter warning

| Warning | Fix |
|---|---|
| `Missing docstring for Symbol` | Add a Julia docstring above the definition — use `/docstrings` |
| `@docs block has unknown object Symbol` | Check spelling and module qualification |
| `Duplicate binding` | Two `@docs` blocks reference the same symbol — remove the duplicate |
| `Cross-reference not found` | Fix the `[text](@ref)` link or add the target docstring |

### Update figures

Figures live in `docs/assets/`. Reference them in markdown as:

```markdown
![Caption](../assets/MyFigure.png)
```

Re-run `docs/make.jl` after adding new images to verify they appear.

## CI

Documentation is built and deployed by `.github/workflows/ci.yml` on pushes to `main`. The deploy target is configured in `docs/make.jl`:

```julia
deploydocs(; repo="https://github.com/albert-de-montserrat/RheologyCalculator.jl", devbranch="main")
```

Check the CI badge in `README.md` for build status.

## Checklist before pushing doc changes

- [ ] `julia --project=docs docs/make.jl` completes with no errors
- [ ] Every exported symbol has a docstring (`checkdocs = :exports` in `make.jl` will error otherwise)
- [ ] New pages are listed in `make.jl`'s `pages` vector
- [ ] New symbols are in the correct `@docs` block in `api.md`
- [ ] Figures are committed to `docs/assets/`
