# Defining composite rheologies

In `RheologyCalculator.jl` we can build composite rheologies with arbitrary configurations in [series](https://en.wikipedia.org/wiki/Maxwell_model), [parallel](https://en.wikipedia.org/wiki/Maxwell_model), or hybrid, with any degree of nesting.

## Material with a configuration in series

Example of a Maxwell visco-elastic model, with a viscous damper of viscosity $\eta=10^20$ and an elastic spring with $G=10$ GPa and $K=46.67$ GPa($nu=0.4$). First we define the invidivual components

```julia
viscous_damper  = LinearViscosity(1e20)
elastic_spring  = Elasticity(10e9, 46.67e9)
```

and then we stitch them together in a `SeriesModel` object:

```julia
maxwell_material = SeriesModel(viscous_damper, elastic_spring)
```

Note that we can put any arbitrary number of individual rheology models together within a single `SeriesModel` object`

## Material with a configuration in parallel

Now let's define the Kelvin-Voigt visco-elastic model, with the same individual elements as in the example above. In this case, we just need to use put them together in a `ParallelModel` object:

```julia
KelvinVoigt_material = ParallelModel(viscous_damper, elastic_spring)
```

As before, we can have any number of rheology models within a single `ParallelModel`

## Hybrid series/parallel material

We can also arbitrarially combine any number of series and parallel rheologies. In the following example we define the Kelvin representation of the [Burgers viscoelastic material](https://en.wikipedia.org/wiki/Burgers_material)

```julia
parallel_element = ParallelModel(viscous_damper, elastic_spring)
Burgers_material = SeriesModel(viscous_damper, elastic_spring, parallel_element)
```

