---
title: 'RheologyCalculator.jl: Composable local rheological models for geodynamics'
tags:
  - Julia
  - geodynamics
  - rheology
  - constitutive models
  - viscoelastoplasticity
  - automatic differentiation
authors:
  - name: Albert de Montserrat
    orcid: 0000-0000-0000-0000 # TODO: confirm ORCID
    affiliation: 1
  - name: Boris J. P. Kaus
    orcid: 0000-0000-0000-0000 # TODO: confirm ORCID
    affiliation: 2
  - name: Thibault Duretz
    orcid: 0000-0000-0000-0000 # TODO: confirm ORCID
    affiliation: 3
affiliations:
  - name: Institute of Geophysics, ETH Zürich, Switzerland # TODO: confirm affiliation
    index: 1
  - name: Institute of Geosciences, Johannes Gutenberg University Mainz, Germany # TODO: confirm affiliation
    index: 2
  - name: Institut für Geowissenschaften, Goethe University Frankfurt, Germany # TODO: confirm affiliation
    index: 3
date: 1 January 2026 # TODO: set submission date
bibliography: paper.bib
---

# Summary

The mechanical behaviour of rocks over geological time is described by
constitutive laws that combine viscous, elastic, and plastic deformation. In
numerical geodynamics these laws are evaluated locally — at every integration
point of a mesh and at every time-step — to relate stress, strain-rate, and
pressure. `RheologyCalculator.jl` is a Julia [@Bezanson2017] package that builds
and solves such *local* rheological models from small, reusable building blocks.
Individual viscous, elastic, and plastic elements are composed in series
(strain-rates add, Maxwell-like), in parallel (stresses add, Kelvin–Voigt-like),
or in arbitrarily nested hybrid networks. From a composed model the package
automatically generates the corresponding nonlinear residual system and solves
it with a Newton–Raphson scheme, returning the updated deviatoric stress,
strain-rate, and pressure at a point.

The design is type-dispatch driven and uses Julia's `@generated` functions
throughout, so that equation assembly and tuple traversals are resolved at
compile time and the resulting solver is allocation-free and type-stable. The
Jacobian required by the Newton solver is obtained by automatic differentiation
with `ForwardDiff.jl` [@Revels2016], and small fixed-size state vectors are
represented with `StaticArrays.jl` [@StaticArrays], making the solver suitable
for the hot inner loop of a Stokes solver. The computational core is
deliberately independent of any particular material catalogue: concrete laws
(e.g. linear and power-law viscosity, diffusion and dislocation creep,
compressible and incompressible elasticity, Drucker–Prager and cap plasticity,
modified Cam-Clay, and rate-and-state friction) are supplied as example
definitions that extend a small state-function interface, and users can add
their own laws the same way.

# Statement of need

Geodynamic Stokes solvers must evaluate the local rheology at every quadrature
point on every time-step, often for strongly nonlinear
visco-elasto-visco-plastic combinations. Traditionally, each combination —
Maxwell, Kelvin–Voigt, Burgers, visco-elasto-plastic (VEP), VEVP,
pressure-dependent yielding, dilatant plasticity — is hand-derived and
hand-coded as a bespoke local stress update. This is tedious and error-prone,
couples the material model tightly to the solver, and makes it difficult to
experiment with new element combinations or to guarantee consistent, converged
local solutions.

`RheologyCalculator.jl` addresses this by treating the local constitutive update
as a graph of composable elements. Given a composite model, it generates the
system of equations automatically, differentiates it exactly via forward-mode
automatic differentiation, and solves it to a prescribed tolerance. This
decouples the *material catalogue* from the *solver machinery*: a new rheology is
added by specialising two dispatch methods and implementing its state functions,
after which it can be freely combined with any existing element in series or
parallel. The package complements the Julia geodynamics ecosystem — such as
`GeoParams.jl` for material parameters [@GeoParams], and the `JustRelax.jl`
[@JustRelax] and `LaMEM` [@Kaus2016] solvers — by providing a general, reusable
engine for the point-wise nonlinear rheological problem that these solvers
repeatedly need to solve.

# Functionality

A model is assembled from elements with `SeriesModel` and `ParallelModel`, which
may be nested:

```julia
using RheologyCalculator
include("rheologies/RheologyDefinitions.jl")

# Maxwell viscoelastic element (viscous + elastic in series)
c = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))

vars   = (; ε = 1.0e-14, θ = 0.0)          # prescribed kinematic inputs
args   = (; τ = 1.0e3, P = 0.0)            # initial guess for unknowns
others = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))  # history / auxiliary fields

x = initial_guess_x(c, vars, args, others)
x = solve(c, x, vars, others)              # Newton solve
```

More elaborate networks — generalized Maxwell/Kelvin–Voigt bodies, Burgers
models, and elastic elements embedded inside parallel branches — are built by
nesting the same two constructors. The solver applies an elastic strain-rate
correction to account for stress history before iterating, and post-processing
utilities recover the updated elastic stress and pressure from the converged
solution vector.

The package is documented with `Documenter.jl` and ships an extensive set of
runnable examples, several of which are validated against analytical solutions
(for example, the mixed Kelvin–Voigt/Maxwell response) and against published
plasticity models, including VEP with a cap [@Popov2025], the hyperbolic yield
surface of @Abbo1995, the modified Cam-Clay formulation of @Golchin2021, and
rate-and-state friction [@Herrendorfer2018]. Correctness is checked by a test
suite covering the equation graphs, Jacobians, and individual rheologies,
together with `Aqua.jl` quality assurance tests.

# Acknowledgements

We acknowledge contributions from the Julia geodynamics community. <!-- TODO:
add funding sources and grant numbers. -->

# References
