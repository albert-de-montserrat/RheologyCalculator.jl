# normalisation factors for the local x vector

""" 
    xnorm = normalisation_x(c::AbstractCompositeModel, char_τ=1.0, char_ε=1.0)

Initial guess for the local solution vector `x` with given characteristic stress and strainrates values
"""
function normalisation_x(c::AbstractCompositeModel, char_τ = 1.0, char_ε = 1.0)
    eqs = generate_equations(c)
    x0 = normalisation_x(eqs, char_τ, char_ε)
    return SVector{length(x0), eltype(x0)}(x0)
end

@generated function normalisation_x(eqs::NTuple{N, CompositeEquation}, char_τ, char_ε) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _normalize_x_value(eqs[i].fn, char_τ, char_ε)
    end
end

@inline _normalize_x_value(::typeof(compute_stress), char_stress, char_strainrate) = char_stress
@inline _normalize_x_value(::typeof(compute_pressure), char_stress, char_strainrate) = char_stress
@inline _normalize_x_value(::typeof(compute_strain_rate), char_stress, char_strainrate) = char_strainrate
@inline _normalize_x_value(::typeof(compute_volumetric_strain_rate), char_stress, char_strainrate) = char_strainrate
@inline _normalize_x_value(::typeof(compute_lambda), char_stress, char_strainrate) = char_stress
@inline _normalize_x_value(::typeof(compute_lambda_parallel), char_stress, char_strainrate) = char_stress
@inline _normalize_x_value(::typeof(compute_plastic_strain_rate), char_stress, char_strainrate) = char_strainrate
@inline _normalize_x_value(::typeof(compute_volumetric_plastic_strain_rate), char_stress, char_strainrate) = char_strainrate

@inline correct_xnorm(::SVector{N,T}, xnorm) where {N, T} = xnorm
@inline correct_xnorm(::SVector{N,T}, ::Nothing) where {N, T} = @SVector ones(T, N)
