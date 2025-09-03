@inline second_invariant(xx, yy, xy) = √(0.5 * (xx^2 + yy^2) + xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)

effective_strain_rate_correction(c::SeriesModel, ε::NTuple{N}, τ0::NTuple{N}, others) where N = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others)

@generated function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε::NTuple{N}, τ0::NTuple{N}, others) where N
    quote
       Base.@ntuple $N i -> effective_strain_rate_correction(c, ε[i], τ0[i], others)
    end
end

@inline effective_strain_rate_correction(::Val{false}, c::SeriesModel, ε::NTuple{N, T}, τ0::NTuple{N}, others) where {N,T} = zero(T)

@inline effective_strain_rate_correction(c::SeriesModel, ε::Number, τ0, others) = effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)

# @generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0, others) where {N}
#     quote
#         i = 0
#         ε_elastic_cor = zero(eltype(ε))
#         Base.@nexprs $N j -> begin
#             if isa(leafs[j], AbstractElasticity)
#                 i += 1
#                 # compute local effect on strainrate tensor due to old elastic stresses
#                 η = compute_viscosity_series(leafs[j], merge((; ε), others))
#                 ε_elastic_cor += τ0[i] / (2 * η)
#             end
#         end
    
#         return ε_elastic_cor
#     end
# end

@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0, others) where {N}
    quote
        @inline
        i = 0
        ε_elastic_cor = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            r = leafs[j]
            i = update_correction_index(r, i)
            η = compute_viscosity(r, merge((; ε), others))
            ε_elastic_cor += effective_strain_rate_correction(r, ε, τ0, others, i)
        end
    
        return ε_elastic_cor
    end
end

effective_strain_rate_correction(c::AbstractRheology, ε, τ0, others, I) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others, I)

effective_strain_rate_correction(::Val{false}, c::AbstractRheology, ε, τ0, others, I) = 0

function effective_strain_rate_correction(::Val{true}, c::AbstractRheology, ε, τ0, others, I)
    η = compute_viscosity(c, merge((; ε), others))
    correction = τ0[I] / (2 * η)
    return correction
end

@inline update_correction_index(c::AbstractRheology, I) = update_correction_index(iselastic(c), I)
@inline update_correction_index(::Val{false}, I) = I
@inline update_correction_index(::Val{true}, I)  = I + 1

@generated function compute_viscosity_parallel(branches::NTuple{N, AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            # compute local effect on strainrate tensor due to old elastic stresses
            η = compute_viscosity_parallel(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        η_eff = 1 / η_eff
        return η_eff
    end
end

@generated function compute_viscosity_series(branches::NTuple{N, AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            # compute local effect on strainrate tensor due to old elastic stresses
            η = compute_viscosity_series(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        η_eff = 1 / η_eff
        return η_eff
    end
end

# elastic trait
@inline iselastic(r::AbstractCompositeModel) = Val(_iselastic(r))
@inline iselastic(::AbstractElasticity)      = Val(true)
@inline iselastic(::AbstractRheology)        = Val(false)

@inline _iselastic(r::AbstractCompositeModel) = _iselastic(r.leafs) || _iselastic(r.branches)

@generated function _iselastic(r::NTuple{N, AbstractRheology}) where N
    quote
        @inline 
        Base.@nexprs $N i -> _iselastic(r[i]) && return true
        return false
    end
end

@inline _iselastic(::AbstractElasticity) = true
@inline _iselastic(::AbstractRheology)   = false
