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
    :compute_vicosity_series,
    :compute_vicosity_parallel,
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
