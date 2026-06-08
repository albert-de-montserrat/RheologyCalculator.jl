fns_state = (
    :compute_strain_rate,
    :compute_stress,
    :compute_volumetric_strain_rate,
    :compute_pressure,
    :compute_lambda,
    :compute_lambda_parallel,
    :compute_plastic_strain_rate,
    :compute_volumetric_plastic_strain_rate,
    :compute_plastic_stress,
    :compute_viscosity,
    :compute_viscosity_series,
    :compute_viscosity_parallel,
)

# ================================================================================================
# ABSTRACT RHEOLOGY (FALLBACK METHODS)
# ================================================================================================
for fn in fns_state
    @eval @inline $fn(r::AbstractRheology; kwargs...) = 0.0e0
end

# ================================================================================================
# WRAPPER FUNCTIONS (for NamedTuple arguments)
# ================================================================================================
for fn in fns_state
    @eval @inline $fn(r::AbstractRheology, kwargs::NamedTuple) = $fn(r; kwargs...)
end

@doc """
    compute_strain_rate(r; τ, kwargs...)
    compute_strain_rate(r, kwargs::NamedTuple)

Return the deviatoric strain-rate contribution of rheology `r` for the supplied
stress-like arguments. Series composites use this as a global/state equation:
the solver finds the stress unknowns whose summed element strain rate matches
the prescribed input `vars.ε`.

Concrete rheologies usually specialize the keyword method. The `NamedTuple`
method forwards to the keyword method and is used by generated residual code.
""" compute_strain_rate

@doc """
    compute_stress(r; ε, kwargs...)
    compute_stress(r, kwargs::NamedTuple)

Return the deviatoric stress contribution of rheology `r` for the supplied
strain-rate-like arguments. Parallel composites use this as a global/state
equation: the solver finds the local strain-rate unknowns whose summed element
stress matches the parent stress.

Concrete rheologies usually specialize the keyword method. The `NamedTuple`
method forwards to the keyword method and is used by generated residual code.
""" compute_stress

@doc """
    compute_volumetric_strain_rate(r; P, kwargs...)
    compute_volumetric_strain_rate(r, kwargs::NamedTuple)

Return the volumetric strain-rate contribution of rheology `r` for the supplied
pressure-like arguments. Volumetric equations are generated only when a
rheology or composite reports volumetric behavior through `_isvolumetric`.
""" compute_volumetric_strain_rate

@doc """
    compute_pressure(r; θ, kwargs...)
    compute_pressure(r, kwargs::NamedTuple)

Return the pressure contribution of rheology `r` for the supplied volumetric
strain-rate-like arguments. This is the volumetric counterpart of
`compute_stress` for parallel composition.
""" compute_pressure

@doc """
    compute_lambda(r; τ, P, λ, kwargs...)
    compute_lambda(r, kwargs::NamedTuple)

Return the local consistency equation for a plastic multiplier `λ` in series
composition. Plastic rheologies specialize this method to connect the yield
function, flow rule, and multiplier unknown.
""" compute_lambda

@doc """
    compute_lambda_parallel(r; τ_pl, P, λ, kwargs...)
    compute_lambda_parallel(r, kwargs::NamedTuple)

Return the local consistency equation for a plastic multiplier `λ` in parallel
composition, usually using branch-local plastic stress variables such as
`τ_pl`.
""" compute_lambda_parallel

@doc """
    compute_plastic_strain_rate(r; τ_pl, λ, P_pl, kwargs...)
    compute_plastic_strain_rate(r, kwargs::NamedTuple)

Return the deviatoric plastic strain-rate residual or contribution for
plastic rheology `r`. The solver treats the associated unknown as `τ_pl`.
""" compute_plastic_strain_rate

@doc """
    compute_volumetric_plastic_strain_rate(r; τ_pl, P_pl, λ, kwargs...)
    compute_volumetric_plastic_strain_rate(r, kwargs::NamedTuple)

Return the volumetric plastic strain-rate residual or contribution for plastic
rheology `r`. This is used by dilatant or volumetric plastic models.
""" compute_volumetric_plastic_strain_rate

@doc """
    compute_plastic_stress(r; τ_pl, kwargs...)
    compute_plastic_stress(r, kwargs::NamedTuple)

Return the stress carried by the plastic branch variable. Plastic rheologies
can specialize this when branch-local plastic stress differs from `τ_pl`
itself.
""" compute_plastic_stress

@doc """
    compute_viscosity(r; kwargs...)
    compute_viscosity(r, kwargs::NamedTuple)

Return an effective viscosity for rheology `r` under the supplied local state.
This helper is used by elastic strain-rate corrections and by user code that
needs a scalar effective viscosity.
""" compute_viscosity

@doc """
    compute_viscosity_series(r; kwargs...)
    compute_viscosity_series(r, kwargs::NamedTuple)

Return an effective viscosity appropriate for series-style aggregation.
Concrete rheologies may specialize this separately from `compute_viscosity`
when the effective-viscosity estimate depends on composition.
""" compute_viscosity_series

@doc """
    compute_viscosity_parallel(r; kwargs...)
    compute_viscosity_parallel(r, kwargs::NamedTuple)

Return an effective viscosity appropriate for parallel-style aggregation.
Concrete rheologies may specialize this separately from `compute_viscosity`
when the effective-viscosity estimate depends on composition.
""" compute_viscosity_parallel

# NOTE: for user defined new functions, add the template below to the appropriate rheology type file (e.g., RheologyDefinitions.jl)
# compute_variable(r::AbstractRheology, kwargs::NamedTuple) = compute_variable(r; kwargs...)
