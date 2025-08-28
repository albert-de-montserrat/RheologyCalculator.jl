# This implements the mode1/mode2 plasticity model proposed in:
# Popov, A. A., Berlie, N., and Kaus, B. J. P.: A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: 
#   stable two-field mixed formulation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-2469, 2025.
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

# ModCamClay ------------------------------------------------------
"""
    ModCamClay{T} <: AbstractPlasticity

Represents a Modified Cam-Clay model, see de Souza Neto book (p. 404)

# Fields
- `M::T`    : The flattening factor for yield (1: circle, < 1: ellipse) = 0.9
- `N::T`    : The flattening factor for potential (1: circle, < 1: ellipse), if N = M, then it's associated
- `r::T`    : The radius along pressure axis = 2e8
- `β::T`    : The asymmetry for compaction cap (1:symmetric, <1 asymmetric)  = 0.1
- `Pt::T`   : The Tensile strength
- `η_vp::T` : The regularisation viscosity
"""
struct ModCamClay{T} <: AbstractPlasticity
    M::T        # flattening factor for yield (1: circle, < 1: ellipse) = 0.9
    N::T        # flattening factor for potential (1: circle, < 1: ellipse)
    r::T        # radius along pressure axis = 2e8
    β::T        # asymmetry for compaction cap (1:symmetric, <1 asymmetric)  = 0.1
    Pt::T       # Tensile strength
    η_vp::T     # regularisation viscosity
end

function ModCamClay(; M=0.9, N=0.9, r=1e8, β=1.0, Pt=-1e7, η_vp=1e20) 
    return ModCamClay(M, N, r, β, Pt, η_vp)
end

#ModCamClay(args...) = ModCamClay(promote(args...)...)
@inline _isvolumetric(::ModCamClay) = true

@inline series_state_functions(::ModCamClay) = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::ModCamClay) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate

@inline function compute_strain_rate(r::ModCamClay; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::ModCamClay; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl/2* (F > -1e-8)
end
@inline function compute_volumetric_strain_rate(r::ModCamClay; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r::ModCamClay; τ_pl = τ, λ = λ, P_pl = P, θ = 0, kwargs...)
    F    = compute_F(r, τ, P)
    return θ_pl* (F > -1e-8) # perhaps this derivative needs to be hardcoded
end

@inline function compute_lambda(r::ModCamClay; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    return -F* (F > -1e-8)  + λ*r.η_vp + λ*1        # last term is for regularisation below yield
end

function compute_F(r::ModCamClay, τII, P)
    (; M, r, β, Pt) = r

    b = P < Pt + r ? one(β) : β

    F  = 1/b *(P - Pt - r)^2  + (τII / M)^2 - r^2 

    # Note that viscoplastic regularisation is taken into account in the residual function
    return F #*(F>-1e-8) 
end

function compute_Q(r::ModCamClay, τII, P) 

    # These parameters are required to compute the constant in the plastic flow
    # potential. Note that this constant does not matter apart when plotting,
    # as we only need derivates of Q in general 
    (; N, r, β, Pt) = r

    b = P < Pt + r ? one(β) : β

    Q  = 1/b *(P - Pt - r)^2  + (τII / N)^2 - r^2 

    return Q
end 

@inline compute_stress(r::ModCamClay; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::ModCamClay; P_pl = 0, kwargs...) = P_pl

@inline function compute_plastic_strain_rate(r::ModCamClay; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ*ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - 0*ε # perhaps this derivative needs to be hardcoded 
end

@inline function compute_volumetric_plastic_strain_rate(r::ModCamClay; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - 0*θ # perhaps this derivative needs to be hardcoded
end
