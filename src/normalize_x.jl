# normalisation factors for the local x vector

""" 
    xnorm = normalisation_x(c::AbstractCompositeModel, char_τ=1.0, char_ε=1.0)

Initial guess for the local solution vector `x` with given characteristic stress and strainrates values
"""
function normalisation_x(c::AbstractCompositeModel, char_τ::T = 1.0, char_ε::T = 1.0) where T
    eqs = generate_equations(c)
    x0 = normalisation_x(eqs, char_τ, char_ε)
    return SA[x0...]
end

@generated function normalisation_x(eqs::NTuple{N, CompositeEquation}, char_τ, char_ε) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _normalize_x_value(eqs[i].fn, char_τ, char_ε)
    end
end

for fn in (:compute_stress, :compute_pressure, :compute_lambda, :compute_lambda_parallel)
    @eval _normalize_x_value(::typeof($fn), char_stress, char_strainrate) = char_stress 
end

for fn in (:compute_strain_rate, :compute_volumetric_strain_rate, :compute_plastic_strain_rate, :compute_volumetric_plastic_strain_rate)
    @eval _normalize_x_value(::typeof($fn), char_stress, char_strainrate) = char_strainrate
end

@inline correct_xnorm(::SVector{N,T}, xnorm) where {N, T} = (@show T; xnorm)
@inline correct_xnorm(::SVector{N,T}, ::Nothing) where {N, T} = @SVector ones(T, N)
