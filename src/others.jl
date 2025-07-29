# complete flatten a tuple
@inline superflatten(t::NTuple{N, Any}) where {N} = superflatten(first(t))..., superflatten(Base.tail(t))...
@inline superflatten(::Tuple{}) = ()
@inline superflatten(x) = (x,)

@generated function _isvolumetric(r::NTuple{N, AbstractRheology}) where {N}
    return quote
        @inline
        b = false
        Base.@nexprs $N i -> b = b || _isvolumetric(r[i])
    end
end

@inline _isvolumetric(::AbstractRheology) = false
@inline _isvolumetric(::Elasticity) = true
@inline _isvolumetric(::BulkElasticity) = true
@inline _isvolumetric(::BulkViscosity) = true
# @inline _isvolumetric(c::AbstractCompositeModel) = _isvolumetric(c.leafs)
@inline _isvolumetric(::Tuple{}) = false

function _isvolumetric(c::AbstractCompositeModel)
    b1 = _isvolumetric(c.leafs)
    branches = c.branches
    b2 = ntuple(Val(length(branches))) do i
        @inline
        _isvolumetric(branches[i])
    end
    b = (b1, b2) |> superflatten
    return reduce(|, b)
end

isvolumetric(c::AbstractCompositeModel) = Val(_isvolumetric(c))

# @generated function harmonic_average(r::NTuple{N, AbstractRheology}, fn::F, args) where {N, F}
#     quote
#         v = 0e0
#         Base.@ntuple $N i -> v += begin
#             x = inv( fn(r[i], args) )
#             x = isinf(x) ? 0e0 : x
#         end
#         return inv(v)
#     end
# end

# harmonic_average_stress(r, args) = harmonic_average(r, compute_stress, args)
# harmonic_average_strain_rate(r, args) = harmonic_average(r, compute_strain_rate, args)
