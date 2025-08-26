# ================================================================================================
# ABSTRACT RHEOLOGY (FALLBACK METHODS)
# ================================================================================================

@inline compute_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_stress(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_volumetric_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_pressure(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_lambda(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_lambda_parallel(r::AbstractRheology; kwargs...) = 0.0e0
@inline compute_plastic_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_volumetric_plastic_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_plastic_stress(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_viscosity(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method

# ================================================================================================
# WRAPPER FUNCTIONS (for NamedTuple arguments)
# ================================================================================================

# splatter wrappers
@inline compute_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_strain_rate(r; kwargs...)
@inline compute_stress(r::AbstractRheology, kwargs::NamedTuple) = compute_stress(r; kwargs...)
@inline compute_volumetric_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_volumetric_strain_rate(r; kwargs...)
@inline compute_pressure(r::AbstractRheology, kwargs::NamedTuple) = compute_pressure(r; kwargs...)
@inline compute_lambda(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda(r; kwargs...)
@inline compute_lambda_parallel(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda_parallel(r; kwargs...)
@inline compute_plastic_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_plastic_strain_rate(r; kwargs...)
@inline compute_volumetric_plastic_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_volumetric_plastic_strain_rate(r; kwargs...)
@inline compute_plastic_stress(r::AbstractRheology, kwargs::NamedTuple) = compute_plastic_stress(r; kwargs...)
@inline compute_viscosity(r::AbstractRheology, kwargs::NamedTuple) = compute_viscosity(r; kwargs...)
