abstract type AbstractCompositeModel  end

# include("rheology_types.jl")
# include("state_functions.jl")
# include("kwargs.jl")
# include("others.jl")
# include("composite.jl")

################ BASIC STRUCTS

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

Base.show(io::IO, ::SeriesModel) = print(io, "SeriesModel")

struct ParallelModel{L, B} <: AbstractCompositeModel # not 100% about the subtyping here, lets see
    leafs::L     # horizontal stacking
    branches::B  # vertical stacking

    function ParallelModel(c::Vararg{Any, N}) where {N}
        leafs = parallel_leafs(c)
        branches = parallel_branches(c)
        return new{typeof(leafs), typeof(branches)}(leafs, branches)
    end
end

Base.show(io::IO, ::ParallelModel) = print(io, "ParallelModel")

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
