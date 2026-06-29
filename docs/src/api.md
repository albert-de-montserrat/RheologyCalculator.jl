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
RheologyCalculator.solve_timed
effective_strain_rate_correction
RheologyCalculator.second_invariant
```

## State Function Interface

Concrete rheology elements participate in generated residual equations by
extending these methods.

```@docs
RheologyCalculator.compute_strain_rate
RheologyCalculator.compute_stress
RheologyCalculator.compute_volumetric_strain_rate
RheologyCalculator.compute_pressure
RheologyCalculator.compute_lambda
RheologyCalculator.compute_lambda_parallel
RheologyCalculator.compute_plastic_strain_rate
RheologyCalculator.compute_volumetric_plastic_strain_rate
RheologyCalculator.compute_plastic_stress
RheologyCalculator.compute_viscosity
RheologyCalculator.compute_viscosity_series
RheologyCalculator.compute_viscosity_parallel
```

## Argument Helpers

```@docs
RheologyCalculator.history_kwargs
RheologyCalculator.differentiable_kwargs
RheologyCalculator.residual_kwargs
RheologyCalculator.generate_args_template
RheologyCalculator.extract_local_kwargs
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
RheologyCalculator.evaluate_state_function
RheologyCalculator.superflatten
RheologyCalculator.isvolumetric
RheologyCalculator.iselastic
```
