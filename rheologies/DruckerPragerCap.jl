# This implements the mode1/mode2 plasticity model proposed in:
# Popov, A. A., Berlie, N., and Kaus, B. J. P.: A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: 
#   stable two-field mixed formulation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-2469, 2025.
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

# DruckerPragerCap ------------------------------------------------------
"""
    DruckerPragerCap{T} <: AbstractPlasticity

Represents a Drucker-Prager plasticity model with cap for pressure-dependent yielding of mode-1 and mode-2 plasticity, 
as described in Popov et al. (2025), Geoscientific Model Development

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
- `η_vp::T`: The Duvaut-Lions regeularisation viscosity for the plasticity model.
- `Pt::T`: The tensile strength (should be < 0).
"""
struct DruckerPragerCap{T} <: AbstractPlasticity
    C::T
    ϕ::T        # in degrees for now
    ψ::T        # in degrees for now
    η_vp::T     # regularisation viscosity
    Pt::T       # Tensile strength

    # computational parameters (precomputed, to speed up later calculations)
    sinϕ::T     # Friction angle
    cosϕ::T     # Friction angle
    sinΨ::T     # Dilation angle
    cosΨ::T     # Dilation angle

    k ::T
    kq::T
    c ::T 
    a ::T
    b ::T
    py::T 
    Ry::T 
    pd::T 
    τd::T 
    pq::T 
end

function DruckerPragerCap(; C=10e6, ϕ=30.0, ψ=0.0, η_vp=1e20, Pt=-1e5) 
    sinϕ = sind(ϕ) # Friction angle
    cosϕ = cosd(ϕ) # Friction angle
    sinΨ = sind(ψ) # Dilation angle
    cosΨ = cosd(ψ) # Dilation angle
    k  = sinϕ
    kq = sinΨ
    c  = C*cosϕ
    a  = sqrt(1.0 + k^2)
    b  = sqrt(1.0 + kq^2)
    py = (Pt + c/a)/(1-k/a)
    Ry = py - Pt
    pd = py - Ry*k/a
    τd = k*pd + c
    pq = pd + kq*τd
    return DruckerPragerCap(C, ϕ, ψ, η_vp, Pt, sinϕ, cosϕ, sinΨ, cosΨ, k, kq, c, a, b, py, Ry, pd, τd, pq)
end

#DruckerPragerCap(args...) = DruckerPragerCap(promote(args...)...)
@inline _isvolumetric(::DruckerPragerCap) = true

@inline series_state_functions(::DruckerPragerCap) = (compute_strain_rate, compute_lambda, compute_volumetric_strain_rate)
@inline parallel_state_functions(::DruckerPragerCap) = compute_stress, compute_pressure, compute_lambda, compute_plastic_strain_rate, compute_volumetric_plastic_strain_rate

@inline function compute_strain_rate(r::DruckerPragerCap; τ = 0, λ = 0, P = 0, kwargs...)
    ε_pl = compute_plastic_strain_rate(r::DruckerPragerCap; τ_pl = τ, λ = λ, P_pl = P, kwargs...)
    F = compute_F(r, τ, P)
    return ε_pl/2* (F > -1e-8)
end
@inline function compute_volumetric_strain_rate(r::DruckerPragerCap; τ = 0, λ = 0, P = 0, kwargs...)
    θ_pl = compute_volumetric_plastic_strain_rate(r::DruckerPragerCap; τ_pl = τ, λ = λ, P_pl = P, θ = 0, kwargs...)
    F    = compute_F(r, τ, P)
    return θ_pl* (F > -1e-8) # perhaps this derivative needs to be hardcoded
end

@inline function compute_lambda(r::DruckerPragerCap; τ = 0, λ = 0, P = 0, kwargs...)
    F = compute_F(r, τ, P)
    return -F* (F > -1e-8)  + λ*r.η_vp + λ*1        # last term is for regularisation below yield
end

# special plastic helper functions
function ismode2_yield(v::DruckerPragerCap{_T}, τII::_T1, P::_T2)  where {_T,_T1,_T2}
    py, τd, pd = v.py, v.τd, v.pd
    return τII*(py - pd) >= τd*(py - P)
end
function ismode2_flowpotential(v::DruckerPragerCap{_T}, τII::_T1, P::_T2)  where {_T,_T1,_T2}
    pq, τd, pd = v.pq, v.τd, v.pd
    return τII*(pq - pd) >= τd*(pq - P)
end

function compute_F(r::DruckerPragerCap, τII, P)
    k, c, py, a, Ry = r.k, r.c, r.py, r.a, r.Ry

    if ismode2_yield(r, τII, P)
        # Mode 2
        F = τII - k * (P)  - c # with fluid pressure (set to zero by default)
    else
        # Mode 1
        Rf   = sqrt(τII^2 + (P - py)^2)
        
        F    = a*(Rf - Ry)  
    end

    # Note that viscoplastic regularisation is taken into account in the residual function
    return F #*(F>-1e-8) 
end

function compute_Q(r::DruckerPragerCap, τ, P) 

    # These parameters are required to compute the constant in the plastic flow
    # potential. Note that this constant does not matter apart when plotting,
    # as we only need derivates of Q in general 
    Rf      = r.pq - r.Pt
    sd      = r.c + r.k*r.pd
    normvRf = sqrt((r.pd - r.pq)^2 + sd^2)/Rf
    pdf     = (r.pd - r.pq)/normvRf + r.pq
    sdf     = sd/normvRf

    if ismode2_flowpotential(r, τ, P) 
        cons =  sdf - r.kq*pdf 
        Q    =  τ - r.kq * (P )  - cons
    else 
        cons =  Rf 
        Rq   =  sqrt(τ^2 + (P - r.pq)^2)
        Q    =  r.b*(Rq - cons)  
    end
    return Q
end 

@inline compute_stress(r::DruckerPragerCap; τ_pl = 0, kwargs...) = τ_pl
@inline compute_pressure(r::DruckerPragerCap; P_pl = 0, kwargs...) = P_pl

@inline function compute_plastic_strain_rate(r::DruckerPragerCap; τ_pl = 0, λ = 0, P_pl = 0, ε = 0, kwargs...)
    return λ*ForwardDiff.derivative(x -> compute_Q(r, x, P_pl), τ_pl) - 0*ε # perhaps this derivative needs to be hardcoded 
end

@inline function compute_volumetric_plastic_strain_rate(r::DruckerPragerCap; τ_pl = 0, λ = 0, P_pl = 0, θ = 0, kwargs...)
    return -λ * ForwardDiff.derivative(x -> compute_Q(r, τ_pl, x), P_pl) - 0*θ # perhaps this derivative needs to be hardcoded
end
