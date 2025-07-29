for fun in (:compute_strain_rate, :compute_volumetric_strain_rate)
    @eval @inline _global_series_state_functions(fn::typeof($fun)) = (fn,)
end

@inline _global_series_state_functions(::F) where {F <: Function} = ()

@generated function global_series_state_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _global_series_state_functions(@inbounds(funs[i]))
        Base.IteratorsMD.flatten(f)
    end
end

function global_series_state_functions(c::SeriesModel)
    fns = series_state_functions(c.leafs)
    return global_series_state_functions(fns)
end

@inline series_state_functions(c::NTuple{N, AbstractCompositeModel}) where {N} = series_state_functions(first(c))..., series_state_functions(Base.tail(c))...
@inline series_state_functions(::Tuple{}) = ()
# @inline series_state_functions(c::ParallelModel)                      = flatten_repeated_functions(parallel_state_functions(c.leafs))

@inline global_series_functions(c::AbstractCompositeModel) = series_state_functions(c.leafs) |> flatten_repeated_functions |> global_series_state_functions

function global_series_functions(c::SeriesModel)
    fns_leafs = series_state_functions(c.leafs)
    fns_branches = series_state_functions(c.branches)
    return (fns_leafs..., fns_branches...) |> flatten_repeated_functions |> global_series_state_functions
    # return global_series_state_functions(fns)
end
