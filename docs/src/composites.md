# [Composites](@id composites)

In `RheologyCalculator.jl` we can build composite rheologies with arbitrary
configurations in series, parallel, or hybrid nested networks.

The two main constructors are [`SeriesModel`](@ref) and [`ParallelModel`](@ref).
Both store direct rheology elements as `leafs` and nested composites as
`branches`; this is the structure used internally by [`generate_equations`](@ref).

In a series model, compatible strain-rate contributions are summed and compared
with the prescribed strain-rate input. In a parallel model, compatible stress
contributions are summed and compared with the parent stress. Nested composites
alternate those roles recursively, which is how Maxwell, Kelvin-Voigt, Burgers,
and more general visco-elasto-plastic networks are represented.

## Material with a configuration in series

Example of a Maxwell visco-elastic model, with a viscous damper of viscosity
$\eta=10^{20}$ and an elastic spring with $G=10$ GPa and $K=46.67$ GPa. First
we define the individual components:

```julia
viscous_damper  = LinearViscosity(1e20)
elastic_spring  = Elasticity(10e9, 46.67e9)
```

and then we stitch them together in a `SeriesModel` object:

```julia
maxwell_material = SeriesModel(viscous_damper, elastic_spring)
```

Any number of individual rheology models can be placed inside one
`SeriesModel`.

## Material with a configuration in parallel

Now let's define the Kelvin-Voigt visco-elastic model with the same individual
elements. In this case, put them together in a `ParallelModel`:

```julia
KelvinVoigt_material = ParallelModel(viscous_damper, elastic_spring)
```

As before, a `ParallelModel` can contain any number of rheology elements.

## Hybrid series/parallel material

We can also combine any number of series and parallel rheologies. In the
following example we define the Kelvin representation of the
[Burgers viscoelastic material](https://en.wikipedia.org/wiki/Burgers_material):

```julia
parallel_element = ParallelModel(viscous_damper, elastic_spring)
Burgers_material = SeriesModel(viscous_damper, elastic_spring, parallel_element)
```

The resulting equation layout can be inspected with:

```julia
eqs = generate_equations(Burgers_material)
x_keys(Burgers_material)
```

`generate_equations` returns the generated residual equations in solver-vector
order. `x_keys` returns the corresponding unknown names, which is useful when
choosing an initial `args` tuple or a normalization vector with
[`normalisation_x`](@ref).
