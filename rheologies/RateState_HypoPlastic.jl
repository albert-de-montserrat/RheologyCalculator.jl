# Hypoplastic Rate & State friction formulation
# 
#module RateState_HypoPlastic


using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import RheologyCalculator.isvolumetric

export  RateStateFriction,
        isvolumetric,
        compute_strain_rate, 
        compute_stress,
        update_Ω

"""
    RateStateFriction{T} <: AbstractViscosity

Implements rate and state dependent frictional viscous creeplas. 
Note that this is not a plasticity model, but a viscous model with state and rate dependence. 

# Fields
- `λ::T`    # fluid pressure ratio
- `μ₀::T`   # reference friction coefficient at reference slip velocity V₀
- `V₀::T`   # reference slip velocity
- `a::T`    # rate & state parameter    
- `b::T`    # rate & state parameter    
- `L::T`    # characteristic slip distance      
- `C::T`    # Cohesion
- `D::T`    # Cell size
"""
struct RateStateFriction{T} <: AbstractViscosity
    λ::T    # fluid pressure ratio
    μ₀::T   # reference friction coefficient at reference slip velocity V₀
    V₀::T   # reference slip velocity
    a::T    # rate & state parameter    
    b::T    # rate & state parameter    
    L::T    # characteristic slip distance
    C::T    # Cohesion 
    D::T    # Cell size
end
RateStateFriction(args...) = RateStateFriction(promote(args...)...)
@inline series_state_functions(::RateStateFriction) = (compute_strain_rate,)

# strain rate as a function of stress, state, and pressure  
@inline function compute_strain_rate(r::RateStateFriction; τ = 0, Ω_old=0, P=0, dt=0, kwargs...)

    Vp = 2 * r.V₀ * sinh(max((τ - r.C), 0) / (r.a *P * (1 - r.λ))) * exp(-(r.μ₀ + r.b * Ω_old) / r.a)
   
    ε = Vp / (2 * r.D)
    return ε
end

# stress as a function of strain rate, state, and pressure
@inline function compute_stress(r::RateStateFriction; ε = 0, Ω_old=0, P=0, dt = 0, kwargs...) 
    Vp  = 2 * r.D * ε
    Ω   = update_Ω(r; ε=ε, Ω_old=Ω_old, dt=dt)
    μd  = r.a * asinh(Vp / (2 * r.V₀) * exp((r.μ₀ + r.b * Ω) / r.a))
    τII = P * (1 - r.λ) * μd + r.C 
    return τII
end

# This updates the state variable Ω based on slip velocity Vp and time step dt, which is needed outside the rheology definition
@inline function update_Ω(r::RateStateFriction; ε = 0, Ω_old=0, dt = 0,  kwargs...)
    Vp  = 2 * r.D * ε
    if (Vp * dt / r.L ≤ 1e-6)
        Ω = log(exp(Ω_old) * (1 - Vp * dt / r.L) + r.V₀ * dt / r.L)
    else
        Ω = log(r.V₀ / Vp + (exp(Ω_old) - r.V₀ / Vp) * exp(-Vp * dt / r.L))
    end
    return Ω
end
    
# --------------------------------------------------------------------
