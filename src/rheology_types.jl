abstract type AbstractRheology end
abstract type AbstractPlasticity <: AbstractRheology end # in case we need specialization at some point
abstract type AbstractElasticity <: AbstractRheology end # in case we need specialization at some point
abstract type AbstractViscosity <: AbstractRheology end # in case we need specialization at some point


## METHODS FOR SERIES MODELS - note that functions for these need to be specified in the rheology element definitions
@inline length_state_functions(r::AbstractRheology) = length(series_state_functions(r))
@inline length_state_functions(r::NTuple{N, AbstractRheology}) where {N} = length_state_functions(first(r))..., length_state_functions(Base.tail(r))...
@inline length_state_functions(r::Tuple{}) = ()
@inline series_state_functions(::AbstractRheology) = error("Rheology not defined")

# returns the flattened statefunctions along with NTuples with global & local element numbers
function series_state_functions(r::NTuple{N, AbstractRheology}, num::MVector{N, Int}) where {N}
    statefuns = (series_state_functions(first(r))..., series_state_functions(Base.tail(r))...)
    len = ntuple(i -> length(series_state_functions(r[i])), N)
    statenum = ntuple(i -> val(i, len, num), Val(sum(len)))
    stateelements = ntuple(i -> val_element(i, len), Val(sum(len)))

    return statefuns, statenum, stateelements
end

# does not allocate:
@inline series_state_functions(r::NTuple{N, AbstractRheology}) where {N} = series_state_functions(first(r))..., series_state_functions(Base.tail(r))...
# @inline series_state_functions(::Tuple{}) = ()

## METHODS FOR PARALLEL MODELS
# table of methods needed per rheology
@inline parallel_state_functions(::AbstractRheology) = error("Rheology not defined")

# handle tuples
@inline parallel_state_functions(r::NTuple{N, AbstractRheology}) where {N} = parallel_state_functions(first(r))..., parallel_state_functions(Base.tail(r))...
@inline parallel_state_functions(::Tuple{}) = ()

function parallel_state_functions(r::NTuple{N, AbstractRheology}, num::MVector{N, Int}) where {N}
    statefuns = (parallel_state_functions(first(r))..., parallel_state_functions(Base.tail(r))...)

    len = ntuple(i -> length(parallel_state_functions(r[i])), N)
    statenum = ntuple(i -> val(i, len, num), Val(sum(len)))
    stateelements = ntuple(i -> val_element(i, len), Val(sum(len)))

    return statefuns, statenum, stateelements
end

@generated function flatten_repeated_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> i == 1 ? (funs[1],) : (funs[i] ∉ funs[1:(i - 1)] ? (funs[i],) : ())
        Base.IteratorsMD.flatten(f)
    end
end

function get_unique_state_functions(composite::NTuple{N, AbstractRheology}, model::Symbol) where {N}
    funs = if model === :series
        get_unique_state_functions(composite, series_state_functions)
    elseif model === :parallel
        get_unique_state_functions(composite, parallel_state_functions)
    else
        error("Model not defined. Accepted models are :series or :parallel")
    end
    return funs
end

function get_unique_state_functions(composite::NTuple{N, AbstractRheology}, state_fn) where {N}
    funs = state_fn(composite)
    # get unique state functions
    return flatten_repeated_functions(funs)
end
