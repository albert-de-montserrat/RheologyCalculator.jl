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

```julia-repl
julia> p = ParallelModel(damper2, spring1)
|--⟦▪̲̅▫̲̅▫̲̅▫̲̅¹--|
|--/\/\/¹--|

julia> c = SeriesModel(damper3, spring2, p)
--⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--|--⟦▪̲̅▫̲̅▫̲̅▫̲̅²--|
                    |--/\/\/²--|
```

In this case, the inout variables are the strain rate `ε` and the volumetric strain rate `θ`, which are constant over the time
```julia
vars = (; ε = 1.0e-15, θ = 1.0e-20)  # input variables (constant)
```
Initial guess for the variables we want to solve for:
```julia
args = (; τ = 1.0e3, P = 1.0e6) 
```
Solution vector `x` contains the initial guess for the variables we want to solve for:
```julia
x   = initial_guess_x(c, vars, args, others)
```

```julia
ntime = 200
dt    = 1.0e9
τ1    = zeros(ntime)
τ2    = zeros(ntime)
P1    = zeros(ntime)
P2    = zeros(ntime)
t_v   = zeros(ntime)
τ_e   = (0.0, 0.0)
P_e   = (0.0, 0.0)
t     = 0.0
for i in 2:ntime
    # non-differentiable variables needed to evaluate the state functions
    others = (; dt = dt, τ0 = τ_e, P0 = P_e) 
    # solve the system of equations
    x      = solve(c, x, vars, others)
    # Post-process the results
    τ_e    = compute_stress_elastic(c, x, others)   # elastic stress
    P_e    = compute_pressure_elastic(c, x, others) # elastic pressure
    # Store the results
    τ1[i], τ2[i] = τ_e
    P1[i], P2[i] = P_e
    t     += others.dt
    t_v[i] = t
end
```