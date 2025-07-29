abstract type AbstractRheology end
abstract type AbstractPlasticity <: AbstractRheology end # in case we need spacilization at some point
abstract type AbstractElasticity <: AbstractRheology end # in case we need spacilization at some point
abstract type AbstractViscosity <: AbstractRheology end # in case we need spacilization at some point

"""
    LinearViscosity{T} <: AbstractViscosity

Represents a linear viscosity model following Newton's law of viscosity.

# Fields
- `η::T`: The dynamic viscosity coefficient.
"""
struct LinearViscosity{T} <: AbstractViscosity
    η::T
end
local_kwargs(r::LinearViscosity) = (τ=0,)

"""
    BulkViscosity{T} <: AbstractViscosity

Represents the bulk viscosity of a material. Bulk viscosity is a material property that characterizes resistance to uniform compression or expansion.

# Fields
- `χ::T`: The value of the bulk viscosity.
"""
struct BulkViscosity{T} <: AbstractViscosity
    χ::T
end

"""
    LinearViscosityStress{T} <: AbstractViscosity

Represents a linear viscosity model that operates on stress rather than strain rate.

# Fields
- `η::T`: The dynamic viscosity coefficient.
"""
struct LinearViscosityStress{T} <: AbstractViscosity
    η::T
end

"""
    PowerLawViscosity{T,I} <: AbstractViscosity

Represents a power-law viscosity model where viscosity depends on strain rate.

# Fields
- `η::T`: The viscosity coefficient.
- `n::I`: The power-law exponent (note: not promoted to floating point by default).
"""
struct PowerLawViscosity{T,I} <: AbstractViscosity
    η::T
    n::I # DO NOT PROMOTE TO FP BY DEFAULT
end

"""
    DiffusionCreep{I,T} <: AbstractViscosity

Represents diffusion creep deformation mechanism in materials.

# Fields
- `n::I`: The stress exponent.
- `r::T`: The water fugacity exponent.
- `p::T`: The grain size exponent.
- `A::T`: The material-specific rheological parameter.
- `E::T`: The activation energy.
- `V::T`: The activation volume.
- `R::T`: The universal gas constant.
"""
struct DiffusionCreep{I,T} <: AbstractViscosity
    n::I
    r::T
    p::T
    A::T
    E::T
    V::T
    R::T
end
DiffusionCreep(args...) = DiffusionCreep(args[1], promote(args[2:end]...)...)

"""
    DislocationCreep{I,T} <: AbstractViscosity

Represents dislocation creep deformation mechanism in materials.

# Fields
- `n::I`: The power-law exponent.
- `r::T`: The exponent of water-fugacity.
- `A::T`: The material specific rheological parameter.
- `E::T`: The activation energy.
- `V::T`: The activation volume.
- `R::T`: The universal gas constant.
"""
struct DislocationCreep{I,T} <: AbstractViscosity
    n::I # power-law exponent
    r::T # exponent of water-fugacity
    A::T # material specific rheological parameter
    E::T # activation energy
    V::T # activation volume
    R::T # universal gas constant
end
DislocationCreep(args...) = DislocationCreep(args[1], promote(args[2:end]...)...)

"""
    LTPViscosity{T} <: AbstractViscosity

Represents a low-temperature plasticity (LTP) viscosity model.

# Fields
- `ε0::T`: The reference strain rate (default: 6.2e-13).
- `Q::T`: The activation energy (default: 76).
- `σb::T`: The brittle strength (default: 1.8 GPa).
- `σr::T`: The reference stress (default: 3.4 GPa).
"""
struct LTPViscosity{T} <: AbstractViscosity
    ε0::T # 6.2e-13
    Q::T  # 76
    σb::T # 1.8 GPa
    σr::T # 3.4 Gpa
end
LTPViscosity(args...) = LTPViscosity(promote(args...)...)

"""
    Elasticity{T} <: AbstractElasticity

Represents elastic deformation with both shear and bulk components.

# Fields
- `G::T`: The shear modulus.
- `K::T`: The bulk modulus.
"""
struct Elasticity{T} <: AbstractElasticity
    G::T
    K::T
end

"""
    IncompressibleElasticity{T} <: AbstractElasticity

Represents incompressible elastic deformation (shear only).

# Fields
- `G::T`: The shear modulus.
"""
struct IncompressibleElasticity{T} <: AbstractElasticity
    G::T
end

"""
    BulkElasticity{T} <: AbstractElasticity

Represents bulk elastic deformation (volumetric compression/expansion only).

# Fields
- `K::T`: The bulk modulus.
"""
struct BulkElasticity{T} <: AbstractElasticity
    K::T
end

"""
    DruckerPrager{T} <: AbstractPlasticity

Represents the Drucker-Prager plasticity model for pressure-dependent yielding.

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
"""
struct DruckerPrager{T} <: AbstractPlasticity
    C::T
    ϕ::T # in degrees for now
    ψ::T # in degrees for now
end

DruckerPrager(args...) = DruckerPrager(promote(args...)...)

## METHODS FOR SERIES MODELS
@inline length_state_functions(r::AbstractRheology) = length(series_state_functions(r))
@inline length_state_functions(r::NTuple{N,AbstractRheology}) where {N} = length_state_functions(first(r))..., length_state_functions(Base.tail(r))...
@inline length_state_functions(r::Tuple{}) = ()

# table of methods needed per rheology
#@inline series_state_functions(::LinearViscosity)          = (compute_strain_rate,)
#@inline series_state_functions(::LTPViscosity)             = (compute_strain_rate,)
#@inline series_state_functions(::PowerLawViscosity)        = (compute_strain_rate,)
#@inline series_state_functions(::IncompressibleElasticity) = (compute_strain_rate, )

types = (:LinearViscosity, :LTPViscosity, :DiffusionCreep, :DislocationCreep, :PowerLawViscosity, :IncompressibleElasticity)
for t in types
    @eval @inline series_state_functions(::($t)) = (compute_strain_rate,)

end
@inline series_state_functions(::LinearViscosityStress) = (compute_stress,)
@inline series_state_functions(::BulkViscosity) = (compute_volumetric_strain_rate,)
@inline series_state_functions(::Elasticity) = compute_strain_rate, compute_volumetric_strain_rate
@inline series_state_functions(::BulkElasticity) = (compute_volumetric_strain_rate,)
# @inline series_state_functions(::DruckerPrager) = compute_strain_rate, compute_volumetric_strain_rate, compute_lambda
@inline series_state_functions(::DruckerPrager) = (compute_strain_rate, compute_lambda)
#@inline series_state_functions(r::Series) = series_state_functions(r.elements)
@inline series_state_functions(::AbstractRheology) = error("Rheology not defined")
# handle tuples


# returns the flattened statefunctions along with NTuples with global & local element numbers
function series_state_functions(r::NTuple{N,AbstractRheology}, num::MVector{N,Int}) where {N}
    statefuns = (series_state_functions(first(r))..., series_state_functions(Base.tail(r))...)
    len = ntuple(i -> length(series_state_functions(r[i])), N)
    statenum = ntuple(i -> val(i, len, num), Val(sum(len)))
    stateelements = ntuple(i -> val_element(i, len), Val(sum(len)))

    return statefuns, statenum, stateelements
end

# does not allocate:
@inline series_state_functions(r::NTuple{N,AbstractRheology}) where {N} = series_state_functions(first(r))..., series_state_functions(Base.tail(r))...
# @inline series_state_functions(::Tuple{}) = ()

## METHODS FOR PARALLEL MODELS
# table of methods needed per rheology
@inline parallel_state_functions(::LinearViscosity) = (compute_stress,)
@inline parallel_state_functions(::PowerLawViscosity) = (compute_stress,)
@inline parallel_state_functions(::Elasticity) = compute_stress, compute_pressure
@inline parallel_state_functions(::BulkElasticity) = (compute_pressure,)
@inline parallel_state_functions(::BulkViscosity) = (compute_pressure,)
@inline parallel_state_functions(::IncompressibleElasticity) = (compute_stress,)
@inline parallel_state_functions(::DruckerPrager) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate
@inline parallel_state_functions(::AbstractRheology) = error("Rheology not defined")

# handle tuples
@inline parallel_state_functions(r::NTuple{N,AbstractRheology}) where {N} = parallel_state_functions(first(r))..., parallel_state_functions(Base.tail(r))...
@inline parallel_state_functions(::Tuple{}) = ()

function parallel_state_functions(r::NTuple{N,AbstractRheology}, num::MVector{N,Int}) where {N}
    statefuns = (parallel_state_functions(first(r))..., parallel_state_functions(Base.tail(r))...)

    len = ntuple(i -> length(parallel_state_functions(r[i])), N)
    statenum = ntuple(i -> val(i, len, num), Val(sum(len)))
    stateelements = ntuple(i -> val_element(i, len), Val(sum(len)))

    return statefuns, statenum, stateelements
end

@generated function flatten_repeated_functions(funs::NTuple{N,Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> i == 1 ? (funs[1],) : (funs[i] ∉ funs[1:(i-1)] ? (funs[i],) : ())
        Base.IteratorsMD.flatten(f)
    end
end

function get_unique_state_functions(composite::NTuple{N,AbstractRheology}, model::Symbol) where {N}
    funs = if model === :series
        get_unique_state_functions(composite, series_state_functions)
    elseif model === :parallel
        get_unique_state_functions(composite, parallel_state_functions)
    else
        error("Model not defined. Accepted models are :series or :parallel")
    end
    return funs
end

function get_unique_state_functions(composite::NTuple{N,AbstractRheology}, state_fn) where {N}
    funs = state_fn(composite)
    # get unique state functions
    return flatten_repeated_functions(funs)
end
