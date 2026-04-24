# API

## Composite Models

```@docs
SeriesModel
ParallelModel
CompositeModel
generate_equations
x_keys
```

## Solver

```@docs
initial_guess_x
normalisation_x
solve
effective_strain_rate_correction
RheologyCalculator.second_invariant
```

## Post-Processing

```@docs
RheologyCalculator.compute_stress_elastic
RheologyCalculator.compute_pressure_elastic
```

## Internals

These functions are mostly useful when extending the package or debugging the
generated residual system.

```@docs
RheologyCalculator.CompositeEquation
RheologyCalculator.compute_residual
RheologyCalculator.generate_args_template
RheologyCalculator.extract_local_kwargs
RheologyCalculator.evaluate_state_function
RheologyCalculator.superflatten
RheologyCalculator.isvolumetric
RheologyCalculator.iselastic
RheologyCalculator.history_kwargs
RheologyCalculator.differentiable_kwargs
RheologyCalculator.residual_kwargs
```
