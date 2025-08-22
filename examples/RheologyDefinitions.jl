# Here we define individual rheological elements, which is is not part of the computational core of 
# RheologyCalculator as it may depend on your local implementation.
using RheologyCalculator, ForwardDiff

# These functions need to be imported, as we use multiple dispatch to extend them here
import RheologyCalculator: series_state_functions, parallel_state_functions
import RheologyCalculator: compute_strain_rate, compute_stress, compute_pressure, compute_volumetric_strain_rate, compute_volumetric_plastic_strain_rate
import RheologyCalculator: compute_plastic_strain_rate, compute_plastic_stress, compute_lambda, compute_lambda_parallel, _isvolumetric
import RheologyCalculator: _isvolumetric

# Linear Viscosity ---------------------------------------------------
"""
    LinearViscosity{T} <: AbstractViscosity

Represents a linear viscosity model following Newton's law of viscosity.

# Fields
- `η::T`: The dynamic viscosity coefficient.
"""
struct LinearViscosity{T} <: AbstractViscosity
    η::T
end
@inline series_state_functions(::LinearViscosity) = (compute_strain_rate,)
@inline parallel_state_functions(::LinearViscosity) = (compute_stress,)

@inline compute_strain_rate(r::LinearViscosity; τ = 0, kwargs...) = τ / (2 * r.η)
@inline compute_stress(r::LinearViscosity; ε = 0, kwargs...) = ε * 2 * r.η
# --------------------------------------------------------------------

# Linear Viscosity but only defined as a function of strainrate ------
"""
    LinearViscosityStress{T} <: AbstractViscosity

Represents a linear viscosity model that operates on stress rather than strain rate.

# Fields
- `η::T`: The dynamic viscosity coefficient.
"""
struct LinearViscosityStress{T} <: AbstractViscosity
    η::T
end
@inline compute_stress(r::LinearViscosityStress; ε = 0, kwargs...) = ε * 2 * r.η
# --------------------------------------------------------------------

# PowerLawViscosity --------------------------------------------------
"""
    PowerLawViscosity{T,I} <: AbstractViscosity

Represents a power-law viscosity model where viscosity depends on strain rate.

# Fields
- `η::T`: The viscosity coefficient.
- `n::I`: The power-law exponent (note: not promoted to floating point by default).
"""
struct PowerLawViscosity{T, I} <: AbstractViscosity
    η::T
    n::I # DO NOT PROMOTE TO FP BY DEFAULT
end
@inline series_state_functions(::PowerLawViscosity) = (compute_strain_rate,)
@inline parallel_state_functions(::PowerLawViscosity) = (compute_stress,)

@inline compute_strain_rate(r::PowerLawViscosity; τ = 0, kwargs...) = τ^r.n / (2 * r.η)
@inline compute_stress(r::PowerLawViscosity; ε = 0, kwargs...) = ε^(1 / r.n) * (2 * r.η)^(1 / r.n)
# --------------------------------------------------------------------

# Elasticity ---------------------------------------------------------
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
@inline _isvolumetric(::Elasticity) = true
@inline series_state_functions(::Elasticity) = (compute_strain_rate, compute_volumetric_strain_rate)
@inline parallel_state_functions(::Elasticity) = (compute_stress, compute_pressure)


@inline compute_strain_rate(r::Elasticity; τ = 0, τ0 = 0, dt = 0, kwargs...) = (τ - τ0) / (2 * r.G * dt)
@inline compute_stress(r::Elasticity; ε = 0, τ0 = 0, dt = 0, kwargs...) = τ0 + 2 * r.G * dt * ε
@inline compute_volumetric_strain_rate(r::Elasticity; P = 0, P0 = 0, dt = 0, kwargs...) = (P - P0) / (r.K * dt)
@inline compute_pressure(r::Elasticity; θ = 0, P0 = 0, dt = 0, kwargs...) = P0 + r.K * dt * θ
# --------------------------------------------------------------------

# Bulk Elasticity ----------------------------------------------------
"""
    BulkElasticity{T} <: AbstractElasticity

Represents bulk elastic deformation (volumetric compression/expansion only).

# Fields
- `K::T`: The bulk modulus.
"""
struct BulkElasticity{T} <: AbstractElasticity
    K::T
end
@inline _isvolumetric(::BulkElasticity) = true
@inline series_state_functions(::BulkElasticity) = (compute_volumetric_strain_rate,)
@inline parallel_state_functions(::BulkElasticity) = (compute_pressure,)

@inline compute_volumetric_strain_rate(r::BulkElasticity; P = 0, P0 = 0, dt = 0, kwargs...) = (P - P0) / (r.K * dt)
@inline compute_pressure(r::BulkElasticity; θ = 0, P0 = 0, dt = 0, kwargs...) = P0 + r.K * dt * θ
# --------------------------------------------------------------------

# Bulk Viscosity -----------------------------------------------------
"""
    BulkViscosity{T} <: AbstractViscosity

Represents the bulk viscosity of a material. Bulk viscosity is a material property that characterizes resistance to uniform compression or expansion.

# Fields
- `χ::T`: The value of the bulk viscosity.
"""
struct BulkViscosity{T} <: AbstractViscosity
    χ::T
end
@inline _isvolumetric(::BulkViscosity) = true
@inline series_state_functions(::BulkViscosity) = (compute_volumetric_strain_rate,)
@inline parallel_state_functions(::BulkViscosity) = (compute_pressure,)

@inline compute_volumetric_strain_rate(r::BulkViscosity; P = 0, kwargs...) = P / r.χ
@inline compute_pressure(r::BulkViscosity; θ = 0, kwargs...) = θ * r.χ
# --------------------------------------------------------------------

# IncompressibleElasticity -------------------------------------------
"""
    IncompressibleElasticity{T} <: AbstractElasticity

Represents incompressible elastic deformation (shear only).

# Fields
- `G::T`: The shear modulus.
"""
struct IncompressibleElasticity{T} <: AbstractElasticity
    G::T
end
@inline series_state_functions(::IncompressibleElasticity) = (compute_strain_rate,)
@inline parallel_state_functions(::IncompressibleElasticity) = (compute_stress,)

@inline compute_strain_rate(r::IncompressibleElasticity; τ = 0, τ0 = 0, dt = 0, kwargs...) = (τ - τ0) / (2 * r.G * dt)
@inline compute_stress(r::IncompressibleElasticity; ε = 0, τ0 = 0, dt = 0, kwargs...) = τ0 + 2 * r.G * dt * ε
# --------------------------------------------------------------------

# LTPViscosity -------------------------------------------------------

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
@inline series_state_functions(::LTPViscosity) = (compute_strain_rate,)

@inline compute_strain_rate(r::LTPViscosity; τ = 0, kwargs...) = max(r.ε0 * sinh(r.Q * (τ - r.σb) / r.σr), 0.0)
# @inline compute_strain_rate(r::LTPViscosity; τ = 0, kwargs...) = r.ε0 * sinh(r.Q * (τ - r.σb) / r.σr)
@inline compute_stress(r::LTPViscosity; ε = 0, kwargs...) = r.σr / r.Q * asinh(ε / r.ε0) + r.σb
# --------------------------------------------------------------------

# DruckerPrager ------------------------------------------------------
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

@inline _isvolumetric(::DruckerPrager) = false
# @inline series_state_functions(::DruckerPrager) = (compute_strain_rate, compute_volumetric_strain_rate, compute_lambda)
# @inline parallel_state_functions(::DruckerPrager) = (compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate)
@inline series_state_functions(::DruckerPrager) = (compute_strain_rate, compute_lambda)
@inline parallel_state_functions(::DruckerPrager) = (compute_stress, compute_plastic_strain_rate, compute_lambda_parallel)

@inline function compute_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P_pl = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::DruckerPrager; τ_pl = τ, λ = λ, P_pl = P_pl, kwargs...)
    return ε_pl
end

@inline compute_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl

@inline function compute_volumetric_strain_rate(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    return λ * ForwardDiff.derivative(x -> compute_Q(r, τ, x), P) # perhaps this derivative needs to be hardcoded
end

@inline compute_pressure(r::DruckerPrager; P_pl = 0, kwargs...) = P_pl

@inline function compute_lambda(r::DruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    η_χ = 1.0  # Lagrange multiplier, value doesn't matter
    return F - λ * η_χ * (F < 0)
end

@inline function compute_lambda_parallel(r::DruckerPrager; τ_pl = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ_pl, P)
    η_χ = 1.0  # Lagrange multiplier, value doesn't matter
    return F - λ * η_χ * (F < 0)
end

@inline function compute_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ / 2 * ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) # perhaps this derivative needs to be hardcoded
end

@inline function compute_volumetric_plastic_strain_rate(r::DruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) # perhaps this derivative needs to be hardcoded
end

@inline compute_plastic_stress(r::DruckerPrager; τ_pl = 0, kwargs...) = τ_pl

# Helper functions for DruckerPrager
function compute_F(r::DruckerPrager, τ, P)
    F = (τ - P * sind(r.ϕ) - r.C * cosd(r.ϕ))
    return F * (F > 0)
end
compute_Q(r::DruckerPrager, τ, P) = τ - P * sind(r.ψ)
# --------------------------------------------------------------------

# DiffusionCreep -----------------------------------------------------
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
struct DiffusionCreep{I, T} <: AbstractViscosity
    n::I
    r::T
    p::T
    A::T
    E::T
    V::T
    R::T
end
DiffusionCreep(args...) = DiffusionCreep(args[1], promote(args[2:end]...)...)
@inline series_state_functions(::DiffusionCreep) = (compute_strain_rate,)
@inline parallel_state_functions(::DiffusionCreep) = (compute_stress,)

@inline function compute_strain_rate(c::DiffusionCreep; τ = 0, T = 0, P = 0, f = 0, args...)
    (; n, r, p, A, E, V, R) = c

    # Q = 52
    # ε = A * τ^n * exp(-Q)
    ε = A * τ^n * exp(-(E + P * V) / (R * T))
    return ε
end
# --------------------------------------------------------------------


# DislocationCreep ---------------------------------------------------  
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
struct DislocationCreep{I, T} <: AbstractViscosity
    n::I # power-law exponent
    r::T # exponent of water-fugacity
    A::T # material specific rheological parameter
    E::T # activation energy
    V::T # activation volume
    R::T # universal gas constant
end
DislocationCreep(args...) = DislocationCreep(args[1], promote(args[2:end]...)...)
@inline series_state_functions(::DislocationCreep) = (compute_strain_rate,)
@inline parallel_state_functions(::DislocationCreep) = (compute_stress,)

@inline function compute_strain_rate(c::DislocationCreep; τ = 0, T = 0, P = 0, f = 0, args...)
    (; n, r, A, E, V, R) = c
    # Q = 73
    # ε = A * τ^n * f^r * exp(-Q)
    ε = A * τ^n * f^r * exp(-(E + P * V) / (R * T))
    return ε
end
# --------------------------------------------------------------------


