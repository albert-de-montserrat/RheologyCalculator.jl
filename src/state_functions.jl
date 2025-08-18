# compute_strain_rate methods
@inline compute_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_strain_rate(r; kwargs...) # splatter wrapper

# compute_volumetric_strain_rate methods
@inline compute_volumetric_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
# splatter wrapper
@inline compute_volumetric_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_volumetric_strain_rate(r; kwargs...)

@inline function compute_lambda_parallel(r::DruckerPrager; τ_pl = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ_pl, P)
    η_χ = 1.0  # Lagrange multiplier, value doesn't matter
    return F - λ * η_χ * (F < 0)
end
@inline compute_lambda_parallel(r::AbstractRheology; kwargs...) = 0.0e0
@inline compute_lambda_parallel(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda_parallel(r; kwargs...)

@inline function compute_lambda(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    η_χ = 1.0  # Lagrange multiplier, value doesn't matter
    return F - λ * η_χ * (F < 0)
end

@inline compute_lambda(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
# splatter wrapper
@inline compute_lambda(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda(r; kwargs...)

# special plastic helper functions
function compute_F(r::DruckerPrager, τ, P)
    F = (τ - P * sind(r.ϕ) - r.C * cosd(r.ϕ))
    return F * (F > 0)
end
compute_Q(r::DruckerPrager, τ, P) = τ - P * sind(r.ψ)

# compute_stress methods
@inline compute_stress(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_stress(r::AbstractRheology, kwargs::NamedTuple) = compute_stress(r; kwargs...)  # splatter wrapper

# compute_pressure methods
@inline compute_pressure(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_pressure(r::AbstractRheology, kwargs::NamedTuple) = compute_pressure(r; kwargs...) # splatter wrapper

# compute_plastic_strain_rate methods
@inline function compute_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ / 2 * ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) # perhaps this derivative needs to be hardcoded
end
@inline compute_plastic_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_plastic_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_plastic_strain_rate(r; kwargs...)    # splatter wrapper

# compute_volumetric_plastic_strain_rate methods
@inline function compute_volumetric_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) # perhaps this derivative needs to be hardcoded
end
@inline compute_volumetric_plastic_strain_rate(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_volumetric_plastic_strain_rate(r::AbstractRheology, kwargs::NamedTuple) = compute_volumetric_plastic_strain_rate(r; kwargs...) # splatter wrapper

# compute_plastic_stress method
@inline compute_plastic_stress(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_plastic_stress(r::AbstractRheology, kwargs::NamedTuple) = compute_plastic_stress(r; kwargs...) # splatter wrapper

# plastic consistency condition
@inline compute_lambda(r::AbstractRheology; kwargs...) = 0.0e0 # for any other rheology that doesnt need this method
@inline compute_lambda(r::AbstractRheology, kwargs::NamedTuple) = compute_lambda(r; kwargs...)

