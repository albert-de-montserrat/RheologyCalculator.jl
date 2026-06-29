"""
    superflatten(x)

Recursively flatten nested tuples and return a single flat tuple. Non-tuple
values are wrapped as one-element tuples.
"""
@inline superflatten(t::NTuple{N, Any}) where {N} = superflatten(first(t))..., superflatten(Base.tail(t))...
@inline superflatten(::Tuple{}) = ()
@inline superflatten(x) = (x,)

"""
    isvolumetric(c)

Return `Val(true)` if a rheology or composite contains any volumetric state
functions, otherwise `Val(false)`.
"""
isvolumetric(c::AbstractCompositeModel) = Val(_isvolumetric(c))
isvolumetric(c::AbstractRheology) = Val(_isvolumetric(c))

@generated function _isvolumetric(r::NTuple{N, AbstractRheology}) where {N}
    return quote
        @inline
        b = false
        Base.@nexprs $N i -> b = b || _isvolumetric(r[i])
    end
end

@inline _isvolumetric(::AbstractRheology) = false
@inline _isvolumetric(::Tuple{}) = false

_isvolumetric(c::AbstractCompositeModel) = _isvolumetric(c.leafs, c.branches)

@generated function _isvolumetric(leafs, branches::NTuple{N, Any}) where {N}
    return quote
        @inline
        b1 = _isvolumetric(leafs)
        b2 = Base.@ntuple $N i -> _isvolumetric(branches[i])
        b = (b1, b2) |> superflatten
        return reduce(|, b)
    end
end

