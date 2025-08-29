# using LinearAlgebra
# using StaticArrays
# using ForwardDiff
# using DifferentiationInterface

# abstract type AbstractRheology end
# abstract type AbstractPlasticity <: AbstractRheology end # in case we need spacilization at some point

abstract type AbstractCompositeModel  end

@inline series_state_functions(::AbstractCompositeModel) = ()
@inline parallel_state_functions(::AbstractCompositeModel) = ()

struct CompositeModel{Nstrain, Nstress, T} <: AbstractCompositeModel
    components::T
end

struct SeriesModel{L, B} <: AbstractCompositeModel # not 100% about the subtyping here, lets see
    leafs::L     # horizontal stacking
    branches::B  # vertical stacking

    function SeriesModel(c::Vararg{Any, N}) where {N}
        leafs = series_leafs(c)
        branches = series_branches(c)
        return new{typeof(leafs), typeof(branches)}(leafs, branches)
    end
end


for fun in (:compute_strain_rate, :compute_volumetric_strain_rate)
    @eval @inline _local_series_state_functions(::typeof($fun)) = ()
    @eval @inline _global_series_state_functions(fn::typeof($fun)) = (fn,)
end

@inline _local_series_state_functions(fn::F) where {F <: Function} = (fn,)

@generated function local_series_state_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _local_series_state_functions(@inbounds(funs[i]))
        Base.IteratorsMD.flatten(f)
    end
end

@inline _global_series_state_functions(::F) where {F <: Function} = ()

@generated function global_series_state_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _global_series_state_functions(@inbounds(funs[i]))
        Base.IteratorsMD.flatten(f)
    end
end

struct ParallelModel{L, B} <: AbstractCompositeModel # not 100% about the subtyping here, lets see
    leafs::L     # horizontal stacking
    branches::B  # vertical stacking

    function ParallelModel(c::Vararg{Any, N}) where {N}
        leafs = parallel_leafs(c)
        branches = parallel_branches(c)
        return new{typeof(leafs), typeof(branches)}(leafs, branches)
    end
end

@inline series_leafs(c::NTuple{N, AbstractRheology}) where {N} = c
@inline series_leafs(c::AbstractRheology) = (c,)
@inline series_leafs(::ParallelModel) = ()
@inline series_leafs(::Tuple{}) = ()
@inline series_leafs(c::NTuple{N, Any}) where {N} = series_leafs(first(c))..., series_leafs(Base.tail(c))...

@inline parallel_leafs(c::NTuple{N, AbstractRheology}) where {N} = c
@inline parallel_leafs(c::AbstractRheology) = (c,)
@inline parallel_leafs(::SeriesModel) = ()
@inline parallel_leafs(::Tuple{}) = ()
@inline parallel_leafs(c::NTuple{N, Any}) where {N} = parallel_leafs(first(c))..., parallel_leafs(Base.tail(c))...

@inline series_branches(::NTuple{N, AbstractRheology}) where {N} = ()
@inline series_branches(::AbstractRheology) = ()
@inline series_branches(c::ParallelModel) = (c,)
@inline series_branches(::Tuple{}) = ()
@inline series_branches(c::NTuple{N, Any}) where {N} = series_branches(first(c))..., series_branches(Base.tail(c))...

@inline parallel_branches(::NTuple{N, AbstractRheology}) where {N} = ()
@inline parallel_branches(::AbstractRheology) = ()
@inline parallel_branches(c::SeriesModel) = (c,)
@inline parallel_branches(::Tuple{}) = ()
@inline parallel_branches(c::NTuple{N, Any}) where {N} = parallel_branches(first(c))..., parallel_branches(Base.tail(c))...

Base.size(c::Union{SeriesModel, ParallelModel}) = length(c.leafs), length(c.branches)

for fun in (:compute_stress, :compute_pressure)
    @eval _local_parallel_state_functions(::typeof($fun)) = ()
    @eval @inline _global_parallel_state_functions(fn::typeof($fun)) = (fn,)
end
@inline _local_parallel_state_functions(fn::F) where {F <: Function} = (fn,)

@generated function local_parallel_state_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _local_parallel_state_functions(@inbounds(funs[i]))
        Base.IteratorsMD.flatten(f)
    end
end

@inline _global_parallel_state_functions(::F) where {F <: Function} = ()

@generated function global_parallel_state_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> _global_parallel_state_functions(@inbounds(funs[i]))
        Base.IteratorsMD.flatten(f)
    end
end

@inline series_state_functions(c::NTuple{N, ParallelModel}) where {N} = series_state_functions(first(c))..., series_state_functions(Base.tail(c))...
@inline series_state_functions(c::Tuple{}) = (compute_strain_rate,)

# @inline series_state_functions(c::ParallelModel)                      = flatten_repeated_functions(parallel_state_functions(c.leafs))
@inline series_state_functions(c::ParallelModel) = flatten_repeated_functions(series_state_functions(c.leafs))


# #######################################################################
# # DEAL FIRST WITH THE SERIES PART
# #######################################################################

@inline function global_series_functions(c::SeriesModel)
    fn_leafs = series_state_functions(c.leafs) |> flatten_repeated_functions |> global_series_state_functions
    fn_branches = series_state_functions(c.branches) |> flatten_repeated_functions |> global_series_state_functions
    return (fn_leafs..., fn_branches...) |> flatten_repeated_functions
end
@inline local_series_functions(c::SeriesModel) = series_state_functions(c.leafs) |> flatten_repeated_functions |> local_series_state_functions

@inline global_parallel_functions(c::SeriesModel) = ntuple(i -> parallel_state_functions(c.branches[i].leafs) |> flatten_repeated_functions |> global_parallel_state_functions, Val(count_parallel_elements(c)))
# @inline local_parallel_functions(c::SeriesModel)  = ntuple(i-> parallel_state_functions(c.branches[i].branches) |> flatten_repeated_functions |> local_parallel_state_functions, Val(count_parallel_elements(c)))

@inline function local_parallel_functions(c::SeriesModel)
    Np = count_parallel_elements(c)
    return ntuple(Val(Np)) do i
        branch = c.branches[i]
        Nb = length(branch.branches)
        ntuple(Val(Nb)) do j
            global_series_functions(branch.branches[j])
        end |> Base.IteratorsMD.flatten
    end
end

# simplify working with it
Base.getindex(c::SeriesModel, i::Int) = c.leafs[i]
Base.getindex(c::ParallelModel, i::Int) = c.leafs[i]
