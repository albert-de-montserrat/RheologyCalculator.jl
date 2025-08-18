# compute_strain_rate methods
@inline compute_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_strain_rate(r; kwargs...) # splatter wrapper

# compute_volumetric_strain_rate methods
@inline compute_volumetric_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_volumetric_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_volumetric_strain_rate(r; kwargs...)  # splatter wrapper

# compute_stress methods
@inline compute_stress(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_stress(r::AbstractRheology, kwargs::NamedTuple) = compute_stress(r; kwargs...)  # splatter wrapper

# compute_pressure methods
@inline compute_pressure(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_pressure(r::AbstractRheology, kwargs::NamedTuple) = compute_pressure(r; kwargs...) # splatter wrapper

# compute_plastic_strain_rate methods
@inline compute_plastic_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_plastic_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_plastic_strain_rate(r; kwargs...)    # splatter wrapper

# compute_volumetric_plastic_strain_rate methods
@inline compute_volumetric_plastic_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_volumetric_plastic_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_volumetric_plastic_strain_rate(r; kwargs...) # splatter wrapper

# compute_plastic_stress method
@inline compute_plastic_stress(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_plastic_stress(r::AbstractRheology, kwargs::NamedTuple) = compute_plastic_stress(r; kwargs...) # splatter wrapper

# plastic consistency condition
@inline compute_lambda(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_lambda(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda(r; kwargs...)

# plastic consistency condition
@inline compute_lambda_parallel(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_lambda_parallel(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda(r; kwargs...)
