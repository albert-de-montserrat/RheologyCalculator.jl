# Stress-time curve of a Burger's material

```julia
η1, η2, η3 = 5e19, 1e20, 1e21
damper1    = LinearViscosity(η1)
damper2    = LinearViscosity(η2)
damper3    = LinearViscosity(η3)
G, K1, K2  = 10e9, 46.67e9, 30e9
spring1    = Elasticity(G, K1)
spring2    = Elasticity(G, K2)
```

```julia
p = ParallelModel(damper2, spring1)
c = SeriesModel(damper3, spring2, p)
```

Input variables: strain rate `ε` and volumetric strain rate `θ`, they are constant all the time:
```julia
vars = (; ε = 1.0e-15, θ = 1.0e-20)  # input variables (constant)
```

```julia
args = (; τ = 1.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
```

```julia
others = (; dt = 1.0e10, τ0 = (1.0, 2.0), P0=(0.0, 0.1))       # other non-differentiable variables needed to evaluate the state functions
```

```julia
x   = initial_guess_x(c, vars, args, others)
```