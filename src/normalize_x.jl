# normalisation factors for the local x vector 

""" 
    xnorm = normalisation_x(c::AbstractCompositeModel, char_τ=1.0, char_ε=1.0)

Initial guess for the local solution vector `x` with given characteristic stress and strainrates values
"""
function normalisation_x(c::AbstractCompositeModel,  char_τ=1.0, char_ε=1.0)
    eqs = generate_equations(c)
    x0  = normalisation_x(eqs,  char_τ, char_ε)
    return SVector{length(x0)}(x0)
end

@generated function normalisation_x(eqs::NTuple{N,CompositeEquation}, char_τ, char_ε)  where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _normalize_x_value(eqs[i].fn, char_τ, char_ε )
    end
end

_normalize_x_value(::typeof(compute_stress), char_stress, char_strainrate) = char_stress
_normalize_x_value(::typeof(compute_pressure), char_stress, char_strainrate) = char_stress
_normalize_x_value(::typeof(compute_strain_rate), char_stress, char_strainrate) = char_strainrate
_normalize_x_value(::typeof(compute_volumetric_strain_rate), char_stress, char_strainrate) = char_strainrate
_normalize_x_value(::typeof(compute_lambda), char_stress, char_strainrate) = char_stress
_normalize_x_value(::typeof(compute_lambda_parallel), char_stress, char_strainrate) = char_stress
_normalize_x_value(::typeof(compute_plastic_strain_rate), char_stress, char_strainrate) = char_strainrate



#=
_normalize_x_value(::Val{:λ}, char_stress, char_strainrate) = char_stress
_normalize_x_value(::Val{:τ}, char_stress, char_strainrate) = char_strainrate
_normalize_x_value(::Val{:τ_pl}, char_stress, char_strainrate) = char_strainrate
_normalize_x_value(::Val{:P}, char_stress, char_strainrate) = char_strainrate
_normalize_x_value(::Val{:P_pl}, char_stress, char_strainrate) = char_strainrate
_normalize_x_value(::Val{:ε}, char_stress, char_strainrate) = char_stress
_normalize_x_value(::Val{Any}, char_stress, char_strainrate) = error("Don't know how to normalize value for $(Val{Any})")
=#
