# This implements the mode1/mode2 plasticity model proposed in:
# Popov, A. A., Berlie, N., and Kaus, B. J. P.: A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: 
#   stable two-field mixed formulation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-2469, 2025.
using RheologyCalculator
import ForwardDiff: ForwardDiff
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

# PorousDruckerPrager ------------------------------------------------------
"""
    PorousDruckerPrager{T} <: AbstractPlasticity

Represents a porous Drucker-Prager model

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
- `η_vp::T`: regularisation viscosity (0)
"""
struct PorousDruckerPrager{T} <: AbstractPlasticity
    C::T   # The cohesion parameter.
    ϕ::T   # The friction angle (in degrees).
    ψ::T   # The dilatancy angle (in degrees).
    η_vp   # regularisation viscosity (0)

    # computational parameters (precomputed, to speed up later calculations)
    sinϕ::T     # Friction angle
    cosϕ::T     # Friction angle
    sinΨ::T     # Dilation angle
    cosΨ::T     # Dilation angle
end

function PorousDruckerPrager(; C=5e7, ϕ=35.0, ψ=10.0, η_vp=0.0) 
    sinϕ = sind(ϕ) # Friction angle
    cosϕ = cosd(ϕ) # Friction angle
    sinΨ = sind(ψ) # Dilation angle
    cosΨ = cosd(ψ) # Dilation angle
    return PorousDruckerPrager(C, ϕ, ψ, η_vp, sinϕ, cosϕ, sinΨ, cosΨ)
end

#PorousDruckerPrager(args...) = PorousDruckerPrager(promote(args...)...)
@inline _isvolumetric(::PorousDruckerPrager) = true

@inline series_state_functions(::PorousDruckerPrager) = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::PorousDruckerPrager) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate

@inline function compute_strain_rate(r::PorousDruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::PorousDruckerPrager; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl/2* (F > -1e-8)
end

@inline function compute_volumetric_strain_rate(r::PorousDruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r::PorousDruckerPrager; τ_pl = τ, λ = λ, P_pl = P, θ = 0, kwargs...)
    F    = compute_F(r, τ, P)
    return θ_pl* (F > -1e-8) # perhaps this derivative needs to be hardcoded
end

@inline compute_porosity_rate(r::PorousDruckerPrager; λ = 0, kwargs...) = λ * r.sinψ

@inline function compute_lambda(r::PorousDruckerPrager; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    return -F* (F > -1e-8)  + λ*r.η_vp + λ*1        # last term is for regularisation below yield
end

function compute_F(r::PorousDruckerPrager, τII, p̄, pᶠ)
    (; C, ϕ, ψ, η_vp, sinϕ, cosϕ, sinΨ, cosΨ) = r
  
    f  = τII - C*cosϕ - (p̄ - pᶠ)*sinϕ

    # Note that viscoplastic regularisation is taken into account in the residual function
    return f #*(F>-1e-8) 
end

function compute_Q(r::PorousDruckerPrager, τII, p̄, pᶠ) 

    # These parameters are required to compute the constant in the plastic flow
    # potential. Note that this constant does not matter apart when plotting,
    # as we only need derivates of Q in general 
    (; C, ϕ, ψ, η_vp, sinϕ, cosϕ, sinΨ, cosΨ) = r
     
    q = τII - C*cosϕ - (p̄ - pᶠ)*sinϕ
    
    return q
end 

@inline compute_stress(r::PorousDruckerPrager; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::PorousDruckerPrager; P_pl = 0, kwargs...) = P_pl

@inline function compute_plastic_strain_rate(r::PorousDruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ*ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - 0*ε # perhaps this derivative needs to be hardcoded 
end

@inline function compute_volumetric_plastic_strain_rate(r::PorousDruckerPrager; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - 0*θ # perhaps this derivative needs to be hardcoded
end
