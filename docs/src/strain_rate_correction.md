# Elastic strain-rate correction

When a composite model contains elastic elements, the constitutive equations
include backstress history terms (the stress carried by each spring at the
*previous* time step, ``\boldsymbol{\tau}^o``).  Before passing the prescribed
strain rate ``\boldsymbol{\varepsilon}`` to the Newton-Raphson solver it is
convenient to absorb those terms into an *effective* strain rate
``\boldsymbol{\varepsilon}^{\mathrm{eff}}``, so that the solver always works
with equations of the canonical form

```math
\boldsymbol{\varepsilon}^{\mathrm{eff}} = \frac{\boldsymbol{\tau}}{2\eta^{\mathrm{eff}}}
```

The function [`effective_strain_rate_correction`](@ref) builds
``\boldsymbol{\varepsilon}^{\mathrm{eff}}`` automatically for any composite.
This page explains the algebra behind it, starting from the simplest case and
building up to the fully general treatment.

All elastic elements use a backward-Euler update:

```math
\boldsymbol{\tau}^{\mathrm{el}} = \boldsymbol{\tau}^o + 2G\,\Delta t\,\boldsymbol{\varepsilon}^{\mathrm{el}}
```

so the elastic strain rate is ``(\boldsymbol{\tau}^{\mathrm{el}} - \boldsymbol{\tau}^o)/(2G\Delta t)``.

---

## Maxwell body

**Model:** `SeriesModel(η, Elasticity(G))`.

In a series model strain rates add:

```math
\boldsymbol{\varepsilon} = \frac{\boldsymbol{\tau}}{2\eta}
  + \frac{\boldsymbol{\tau} - \boldsymbol{\tau}^o}{2G\Delta t}
```

Moving the backstress term to the left-hand side so that the right-hand side
depends only on ``\boldsymbol{\tau}``:

```math
\underbrace{\boldsymbol{\varepsilon} + \frac{\boldsymbol{\tau}^o}{2G\Delta t}}
           _{\boldsymbol{\varepsilon}^{\mathrm{eff}}}
= \boldsymbol{\tau}
  \underbrace{\left(\frac{1}{2\eta} + \frac{1}{2G\Delta t}\right)}
             _{1/(2\eta^{\mathrm{eff}})}
```

The Maxwell effective viscosity is the harmonic mean of ``\eta`` and
``G\Delta t``:

```math
\eta^{\mathrm{eff}}_M = \frac{1}{1/\eta + 1/(G\Delta t)}
```

The strain-rate correction is simply ``\boldsymbol{\tau}^o/(2G\Delta t)`` —
the elastic backstress divided by the elastic stiffness.  In code this is
handled by the specialisation of `effective_strain_rate_correction` for
`AbstractElasticity` leafs of the outer `SeriesModel`.

---

## Kelvin-Voigt body

**Model:** `SeriesModel(η₁, ParallelModel(Elasticity(G), η₂))`.

The outer series splits the total strain rate:

```math
\boldsymbol{\varepsilon} = \frac{\boldsymbol{\tau}}{2\eta_1} + \boldsymbol{\varepsilon}^p
```

Inside the parallel (Kelvin-Voigt) branch stresses add:

```math
\boldsymbol{\tau} = 2\eta_2\,\boldsymbol{\varepsilon}^p
  + \left(2G\Delta t\,\boldsymbol{\varepsilon}^p + \boldsymbol{\tau}^o\right)
= (2\eta_2 + 2G\Delta t)\,\boldsymbol{\varepsilon}^p + \boldsymbol{\tau}^o
```

Defining the Kelvin-Voigt effective viscosity ``\eta^{\mathrm{eff}}_{\mathrm{KV}} = \eta_2 + G\Delta t``
and solving for ``\boldsymbol{\varepsilon}^p``:

```math
\boldsymbol{\varepsilon}^p = \frac{\boldsymbol{\tau} - \boldsymbol{\tau}^o}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
```

Substituting back into the series equation and isolating the effective strain
rate:

```math
\boldsymbol{\varepsilon}^{\mathrm{eff}}
= \boldsymbol{\varepsilon} + \frac{\boldsymbol{\tau}^o}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
= \boldsymbol{\tau}\!\left(\frac{1}{2\eta_1} + \frac{1}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}\right)
```

The correction now scales with ``\eta^{\mathrm{eff}}_{\mathrm{KV}}`` rather than
``G\Delta t`` alone, because the elastic backstress must be *distributed across
the stiffness of the whole parallel branch*.  For ``N`` parallel elements the
KV effective viscosity generalises to an arithmetic sum:

```math
\eta^{\mathrm{eff}}_{\mathrm{KV}} = \sum_{i=1}^N \eta_i
```

---

## Mixed Kelvin-Voigt and Maxwell body

**Model:** `SeriesModel(η₁, ParallelModel(η₂, SeriesModel(η₃, Elasticity(G))))`.

The inner `SeriesModel(η₃, Elasticity(G))` is itself a Maxwell element with
its own effective viscosity ``\eta^{\mathrm{eff}}_M = 1/(1/\eta_3 + 1/(G\Delta t))``.
Its strain rate equals the branch strain rate ``\boldsymbol{\varepsilon}^p`` of
the outer parallel block:

```math
\boldsymbol{\varepsilon}^p = \frac{\boldsymbol{\tau}^{\mathrm{Max}}}{2\eta^{\mathrm{eff}}_M}
  - \frac{\boldsymbol{\tau}^o}{2G\Delta t}
\quad\Longrightarrow\quad
\boldsymbol{\tau}^{\mathrm{Max}} = 2\eta^{\mathrm{eff}}_M
  \!\left(\boldsymbol{\varepsilon}^p + \frac{\boldsymbol{\tau}^o}{2G\Delta t}\right)
```

The parallel stress balance is:

```math
\boldsymbol{\tau} = 2\eta_2\,\boldsymbol{\varepsilon}^p + \boldsymbol{\tau}^{\mathrm{Max}}
```

Substituting ``\boldsymbol{\tau}^{\mathrm{Max}}`` and collecting
``\boldsymbol{\varepsilon}^p``:

```math
\boldsymbol{\tau}
= 2(\eta_2 + \eta^{\mathrm{eff}}_M)\,\boldsymbol{\varepsilon}^p
  + \frac{\eta^{\mathrm{eff}}_M}{G\Delta t}\,\boldsymbol{\tau}^o
```

```math
\boldsymbol{\varepsilon}^p
= \frac{1}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
  \!\left(\boldsymbol{\tau} - \eta^{\star}_M\,\boldsymbol{\tau}^o\right),
\qquad
\eta^{\mathrm{eff}}_{\mathrm{KV}} = \eta_2 + \eta^{\mathrm{eff}}_M,
\qquad
\eta^{\star}_M = \frac{\eta^{\mathrm{eff}}_M}{G\Delta t}
```

Substituting back into the outer series equation, the effective strain rate is:

```math
\boldsymbol{\varepsilon}^{\mathrm{eff}}
= \boldsymbol{\varepsilon} + \frac{\eta^{\star}_M\,\boldsymbol{\tau}^o}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
```

The weighting factor ``\eta^{\star}_M \in (0, 1)`` attenuates the backstress:

- **Soft spring** (``G \to 0``): ``\eta^{\star}_M \to 1`` — full backstress, elastic history dominates.
- **Stiff spring** (``G \to \infty``): ``\eta^{\star}_M \to 0`` — backstress vanishes, element behaves as a pure dashpot.

### Example

Under a constant deviatoric strain rate ``\dot{\varepsilon}`` starting from zero stress, this
model obeys a linear first-order ODE whose exact solution is

```math
\tau(t)
= \tau_\infty - \bigl(\tau_\infty - \tau_0\bigr)\,e^{-t/t_{\mathrm{relax}}}
```

with

```math
\tau_0 = \frac{2\eta_1\eta_2}{\eta_1+\eta_2}\,\dot{\varepsilon}
\qquad\text{(initial stress — spring unloaded)}
```

```math
\tau_\infty = \frac{2\eta_1(\eta_2+\eta_3)}{\eta_1+\eta_2+\eta_3}\,\dot{\varepsilon}
\qquad\text{(long-time equilibrium)}
```

```math
t_{\mathrm{relax}} = \frac{(\eta_1+\eta_2)\,\eta_3}{G\,(\eta_1+\eta_2+\eta_3)}
\qquad\text{(relaxation timescale)}
```

At ``t=0`` the spring carries no backstress so only the viscous elements ``\eta_1`` and
``\eta_2`` resist deformation; as ``t \to \infty`` the inner Maxwell dashpot reaches its
steady-state load and ``\tau`` approaches the higher value ``\tau_\infty``.

The full runnable example with convergence plot is in `examples/Maxwell_KV_Maxwell.jl`.
The essential time-stepping structure below follows the same 2D tensor conventions used
throughout the other examples (load `rheologies/RheologyDefinitions.jl` and
`examples/tensor_helpers.jl` before running):

```julia
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

η1 = LinearViscosity(1e22)
η2 = LinearViscosity(1e21)
η3 = LinearViscosity(1e21)
el = IncompressibleElasticity(1e10)

c = SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))

εII    = 1.0e-14
vars   = vars_2D(εII)
args   = (; τ = 2.0e7, P = 0.0)
others = (; dt = 1.0e9, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))

x = initial_guess_x(c, vars, args, others)

τ_char = 2η1.η * (η2.η + η3.η) / (η1.η + η2.η + η3.η) * εII
xnorm  = normalisation_x(c, τ_char, εII)

τ_e = (zero_stress_tensor_2D(),)
P_e = (0.0,)
t   = 0.0
for _ in 1:1_000
    others = (; dt = 1.0e9, τ0 = τ_e, P0 = P_e)
    x      = solve(c, x, vars, others; xnorm0 = xnorm)
    # pass the full x SVector — the elastic element is nested inside SeriesModel(η₃, G)
    # so compute_stress_elastic must look at x[3], not x[1]
    τ_e    = elastic_stress_history_2D(c, x, vars.ε, τ_e, others)
    t     += 1.0e9
end
```

`x_keys(c)` returns `(:τ, :ε, :τ)` for this composite — the three unknowns are the
global stress, the strain rate entering the parallel branch, and the stress inside the
inner Maxwell element. Because the elastic element is nested two levels deep,
`elastic_stress_history_2D` receives the full `x` SVector (not just `x[1]`) so that
`compute_stress_elastic` can locate the spring stress at index 3 and return the correct
backstress for the next step.

---

## Generalized Maxwell body

The two cases above are instances of a general pattern.  Consider any
`SeriesModel` whose branches are `ParallelModel` blocks, each of which may
itself contain Maxwell `SeriesModel` sub-branches.  The effective strain-rate
correction contributed by one `ParallelModel` branch is:

```math
\Delta\boldsymbol{\varepsilon}^{\mathrm{eff}}
= \frac{1}{2\eta_{\mathrm{KV}}} \sum_{i} \eta^{\star}_i\,\boldsymbol{\tau}^o_i
```

where the sum runs over every elastic source (direct elastic leaf or Maxwell
sub-branch) inside the parallel block, and:

```math
\eta_{\mathrm{KV}} = \sum_i \eta^{\mathrm{eff}}_i
  \qquad\text{(arithmetic sum of all branch effective viscosities)}
```

```math
\eta^{\star}_i =
\begin{cases}
  1 & \text{direct elastic leaf of the } \texttt{ParallelModel} \\[4pt]
  \dfrac{\eta^{\mathrm{eff}}_{M_i}}{G_i\Delta t} = \dfrac{\eta_{\nu_i}}{\eta_{\nu_i} + G_i\Delta t}
    & \text{Maxwell } \texttt{SeriesModel} \text{ sub-branch}
\end{cases}
```

For a purely viscous branch ``\eta^{\star} = 1`` but ``\boldsymbol{\tau}^o = 0``
(no elastic history), so its contribution vanishes automatically.

The total effective strain rate seen by the outer solver is the sum of the
existing Maxwell-leaf correction and all branch corrections:

```math
\boldsymbol{\varepsilon}^{\mathrm{eff}}
= \boldsymbol{\varepsilon}
  + \underbrace{\sum_{k \in \text{elastic leafs}} \frac{\boldsymbol{\tau}^o_k}{2G_k\Delta t}}_{\text{Maxwell correction}}
  + \underbrace{\sum_{j \in \text{branches}}
      \frac{1}{2\eta^{(j)}_{\mathrm{KV}}}
      \sum_{i \in \text{branch } j} \eta^{\star}_{i}\,\boldsymbol{\tau}^o_i}_{\text{KV / generalized Maxwell correction}}
```

---

## Implementation

The correction is computed by `effective_strain_rate_correction` (defined in
`src/strain_rate_correction.jl`) and called automatically by `solve` before
each Newton step.  The decomposition into *leafs* and *branches* mirrors the
internal structure of `SeriesModel`:

| Function | Role |
|---|---|
| `effective_strain_rate_correction(leafs, (), ε, τ0, others)` | Maxwell correction for direct elastic leafs |
| `_kv_corrections(branches, ε, τ0, others, offset)` | KV / generalized Maxwell correction for all `ParallelModel` branches |
| `_kv_branch_correction(branch, ε, τ0, others, el_idx_start)` | Correction for one `ParallelModel` branch |
| `_η_KV(leafs, subs, args)` | ``\eta_{\mathrm{KV}} = \sum \eta_i`` for one branch |
| `_η_eff_maxwell(leafs, args)` | ``\eta^{\mathrm{eff}}_M`` (harmonic mean) for a Maxwell sub-branch |
| `_η_eff_elastic(leafs, args)` | ``G\,\Delta t`` for the elastic leaf of a sub-branch |
| `_weighted_backstress(leafs, subs, ε, τ0, args, el_idx_start)` | ``\sum \eta^{\star}_i\,\boldsymbol{\tau}^o_i`` for one branch |

All τ0 index arithmetic and ``\eta^{\star}`` factors are resolved at
*compile time* via `@generated` functions, so the emitted code is a flat
sequence of multiply-adds with no runtime dispatch or branching.

### Worked example

```julia
using RheologyCalculator
include("rheologies/RheologyDefinitions.jl")

η1 = LinearViscosity(1e22)   # outer series dashpot
η2 = LinearViscosity(1e21)   # parallel dashpot (KV branch)
η3 = LinearViscosity(5e20)   # inner Maxwell dashpot
el = IncompressibleElasticity(30e9)

# SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))
c = SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))

dt   = 1e10
τ0   = (1e6,)          # backstress for the one elastic element
ε    = (1e-14,)
others = (; dt, τ0, P0 = (0.0,))

# Analytical expected value
G_dt    = el.G * dt
η_eff_M = 1 / (1/η3.η + 1/G_dt)   # Maxwell eff. viscosity of inner branch
η_star  = η_eff_M / G_dt           # attenuation weight
η_KV    = η2.η + η_eff_M           # KV eff. viscosity of parallel block
correction = η_star * τ0[1] / (2η_KV)

# Computed value
import RheologyCalculator: effective_strain_rate_correction
cor = effective_strain_rate_correction(c, ε, τ0, others)

@assert isapprox(only(cor), correction; rtol = 1e-10)
```
