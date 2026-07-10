# Contributing to RheologyCalculator.jl

Contributions are welcome — bug reports, feature requests, documentation
improvements, new rheologies, and code.

## Reporting issues and getting support

- Please open an [issue](https://github.com/albert-de-montserrat/RheologyCalculator.jl/issues)
  for bugs, questions, or feature requests.
- When reporting a bug, include a minimal reproducible example, the error
  message, and your Julia and package versions (`] status`).

## Contributing code

1. Fork the repository and create a feature branch.
2. Install the development environment:
   ```julia
   using Pkg
   Pkg.develop(path = ".")
   Pkg.instantiate()
   ```
3. Make your change. New rheologies are added by subtyping `AbstractViscosity`,
   `AbstractElasticity`, or `AbstractPlasticity`, specialising
   `series_state_functions` / `parallel_state_functions`, and implementing the
   corresponding state functions (see `rheologies/RheologyDefinitions.jl` for
   canonical examples and `CLAUDE.md` for the architecture).
4. Add or update tests under `test/` and run the suite:
   ```julia
   using Pkg
   Pkg.test("RheologyCalculator")
   ```
5. Open a pull request describing the change. CI must pass before review.

## Code style

Follow the conventions of the surrounding code. The core relies on
`@generated` functions and type dispatch to stay allocation-free — keep new
core code type-stable and avoid runtime allocation in the solver hot path.

## Code of conduct

Please be respectful and constructive in all interactions. By participating you
agree to uphold a welcoming and harassment-free environment.
