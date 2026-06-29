"""
    second_invariant(a)
    second_invariant(xx, yy, xy)
    second_invariant(xx, yy, zz, yz, xz, xy)

Return `a` for a scalar invariant, or the second invariant of a 2D or 3D
symmetric deviatoric tensor stored in Voigt-like component order.
"""
@inline second_invariant(a::Number) = a
@inline second_invariant(xx, yy, xy) = √((xx^2 + yy^2 + (-xx - yy)^2) / 2 +  xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)
@inline second_invariant_value(a::Number) = second_invariant(a)
@inline second_invariant_value(a::NTuple) = second_invariant(a...)

"""
    effective_strain_rate_correction(c, ε, τ0, others)

Return the strain-rate correction due to previous elastic stress history in
composite `c`. Each elastic element contributes `τ0 / (2η)` evaluated at the
current effective viscosity; non-elastic elements contribute zero.

For tensor strain-rate inputs, `ε` and the return value are `NTuple`s in Voigt
order. For a scalar input, the return value is a scalar. [`solve`](@ref) applies
this correction automatically before Newton iteration; call it directly only when
you need the raw correction value.

# Example
```julia
c      = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
ε      = (1e-15, -1e-15, 0.0)            # 2-D Voigt deviatoric tensor
others = (; dt = 1e10, τ0 = (0.0,), P0 = (0.0,))
δε     = effective_strain_rate_correction(c, ε, others.τ0, others)
```
"""
effective_strain_rate_correction(c::SeriesModel, ε::NTuple, τ0::NTuple, others) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others)

# @generated function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε::NTuple, τ0::NTuple{N}, others) where N
#     quote
#         @inline 
#         Base.@ntuple $N i -> effective_strain_rate_correction(c, ε, τ0[i], others)
#     end
# end
function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε::NTuple, τ0::NTuple, others)
    effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)
end

@inline effective_strain_rate_correction(::Val{false}, c::SeriesModel, ε::NTuple{N, T}, τ0::NTuple{N}, others) where {N,T} = zero(T)

@inline effective_strain_rate_correction(c::SeriesModel, ε, τ0, others) = effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)

@inline function effective_strain_rate_correction(leafs::NTuple{N, Any}, branches::NTuple{Nb, Any}, ε, τ0::NTuple{Nτ}, others) where {N, Nb, Nτ}
    return effective_strain_rate_correction(leafs, (), ε, τ0, others)
end

@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0::NTuple{Nτ}, others) where {N, Nτ}
    quote
        @inline
        i = 0
        ε_elastic_cor = ε .* 0
        Base.@nexprs $N j -> begin
            i = update_correction_index(leafs[j], i)
            if i > 0
                # η = compute_viscosity(leafs[j], merge((; ε), others))
                ε_elastic_cor = ε_elastic_cor .+ effective_strain_rate_correction(leafs[j], ε, τ0[i], others, i)
            end
        end
        
        return ε_elastic_cor
    end
end

@inline effective_strain_rate_correction(c::AbstractRheology, ε, τ0, others, I) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others, I)
@inline effective_strain_rate_correction(::Val{false}, c::AbstractRheology, ε, τ0, others, I) = 0

@inline function effective_strain_rate_correction(::Val{true}, c::AbstractRheology, ε, τ0, others, I)
    η = compute_viscosity(c, merge((; ε), others))
    correction = @. τ0 / (2 * η)
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

"""
    iselastic(r)

Return `Val(true)` when `r` is an elastic rheology or a composite containing an
elastic rheology, otherwise `Val(false)`.
"""
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
