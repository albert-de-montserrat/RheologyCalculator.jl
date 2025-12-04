# This implements the mode1/mode2 plasticity model proposed in:
# Popov, A. A., Berlie, N., and Kaus, B. J. P.: A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: 
#   stable two-field mixed formulation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-2469, 2025.
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

# Golchin ------------------------------------------------------
"""
    Golchin{T} <: AbstractPlasticity

Represents a extended Cam-Clay model, see Golchin et al. (2021)

# Fields
- ` M::T`   : "friction":  6*sind(ϕ) / (3 - sind(ϕ)) for tri-axial compression (0.9)
- ` N::T`   : "dilatancy": 6*sind(ψ) / (3 - sind(ψ)) for tri-axial compression (0.6)
- ` γ::T`   : shape factor (0.5, 0: ellipse)
- ` α::T`   : shape factor (0.5, 0: ellipse)
- ` β::T`   : asymmetry for compaction cap (0)
- ` Pt::T`  : Tensile pressure (-1e7)
- ` Pc::T`  : Compaction pressure (3e8)
- ` η_vp::T`: regularisation viscosity (0)
"""
struct Golchin{T} <: AbstractPlasticity
    M::T        # "friction"
    N::T        # "dilatancy"
    γ::T        # shape factor (0: ellipse)
    α::T        # shape factor (0: ellipse)
    β::T        # asymmetry for compaction cap 
    Pt::T       # Tensile pressure
    Pc::T       # Compaction pressure
    η_vp::T     # regularisation viscosity
end

@inline Af(p, pc, pt, γ) = (pc - pt)/(2*π) *(2*atan(γ*(pc+pt-2p)/(2*pc))+π)
@inline Bf(p, pc, pt, M, C, α) = M*C*exp(α*(p - C)/(pc - pt))
@inline Cf(pc, pt, γ) = (pc - pt)/π * atan(γ/2) + (pc + pt)/2  

function Golchin(; M=0.9, N=0.6, γ=0.5, α=0.5, β=0.0, Pt=-1e7, Pc=3e8, η_vp=0.0) 
    return Golchin(M, N, γ, α, β, Pt, Pc, η_vp)
end

#Golchin(args...) = Golchin(promote(args...)...)
@inline _isvolumetric(::Golchin) = true

@inline series_state_functions(::Golchin) = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::Golchin) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate

@inline function compute_strain_rate(r::Golchin; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::Golchin; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl/2* (F > -1e-8)
end
@inline function compute_volumetric_strain_rate(r::Golchin; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r::Golchin; τ_pl = τ, λ = λ, P_pl = P, θ = 0, kwargs...)
    F    = compute_F(r, τ, P)
    return θ_pl* (F > -1e-8) # perhaps this derivative needs to be hardcoded
end

@inline function compute_lambda(r::Golchin; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    return -F* (F > -1e-8)  + λ*r.η_vp + λ*1        # last term is for regularisation below yield
end

function compute_F(r::Golchin, τII, P)
    (; M, N, γ, α, β, Pt, Pc) = r

    C = Cf(Pc, Pt, γ)           
    B = Bf(P, Pc, Pt, M, C, α)
    A = Af(P, Pc, Pt, γ)    

    F  = (P - C)^2/A^2 + (τII - β*P)^2/B^2 - 1

    # Note that viscoplastic regularisation is taken into account in the residual function
    return F #*(F>-1e-8) 
end

function compute_Q(r::Golchin, τII, P) 

    # These parameters are required to compute the constant in the plastic flow
    # potential. Note that this constant does not matter apart when plotting,
    # as we only need derivates of Q in general 
    (; M, N, γ, α, β, Pt, Pc) = r

    C = Cf(Pc, Pt, γ)         
    B = Bf(P, Pc, Pt, N, C, α)
    A = Af(P, Pc, Pt, γ)      

    Q  = (P - C)^2/A^2 + (τII - β*P)^2/B^2 - 1

    return Q
end 

@inline compute_stress(r::Golchin; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::Golchin; P_pl = 0, kwargs...) = P_pl

@inline function compute_plastic_strain_rate(r::Golchin; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ*ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - 0*ε # perhaps this derivative needs to be hardcoded 
end

@inline function compute_volumetric_plastic_strain_rate(r::Golchin; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - 0*θ # perhaps this derivative needs to be hardcoded
end
