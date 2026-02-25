# This gives an example of using a rate and state hypoplasticity model in a 0D setting.
# This is similar to the way rate & state friction is implemented in 
#   Herrendörfer, R., Gerya, T., van Dinther, Y., 2018. An Invariant Rate- and State-Dependent Friction 
#      Formulation for Viscoeastoplastic Earthquake Cycle Simulations. Journal of Geophysical Research: Solid Earth 123, 5018–5051. https://doi.org/10.1029/2017JB015225
#
# It also implements an adaptive timestepping scheme that is crucial

using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie
import Statistics: mean

include("../rheologies/RheologyDefinitions.jl")

# This implements rate & state in an explicit manner, where "state" is not included in the local iterations
include("../rheologies/RateState_HypoPlastic.jl")

"""
    dt = compute_dt_ratestate(rs::RateState_HypoPlastic, re::AbstractElasticity, P, Vp, Ω; dt_min=1e-4, dt_max=1e6, f_courant=1e-3)

Adaptive timestepping function for a rate and state friction model, that takes care of healing and weakening constraints.
"""
function compute_dt_ratestate(rs::RateStateFriction, re::AbstractElasticity, P, Vp, Ω; dt_min=1e-4, dt_max=1e6, f_courant=1e-3)
   
    dt_c = dt_courant(Vp, rs.D; f=f_courant)
    dt_h = dt_healing(rs, Ω)
    θmax = max_state_change(rs, re, P)
    dt_w = dt_weakening(rs, θmax, Vp)
    dt = min(dt_c, dt_h, dt_w) 

    # Upper & lower cutoffs
    dt = clamp(dt, dt_min, dt_max)

    return dt
end


"""
    dt = dt_healing(r::RateStateFriction, state)
Compute time step based on healing criterion
"""
dt_healing(r::RateStateFriction, Ω) =  0.2 * r.L / (r.V₀ * exp(-Ω))

"""
    dt = dt_weakening(r::RateStateFriction, θmax, Vp)
Compute time step based on weakening criterion, `θmax` is maximum allowed state variable change
"""
dt_weakening(r::RateStateFriction, θmax, Vp) =  θmax * r.L / Vp

"""
    dt = dt_courant(Vmax, Δ; f=1e-3)
Compute time step based on Courant criterion, `Vmax` is maximum velocity, `Δ` is the (minimum) spatial resolution, and `f` is a safety factor
"""
dt_courant(Vmax, Δ; f=1e-3) = f * Δ / Vmax


"""
    θmax = max_state_change(rs::RateStateFriction, re::AbstractElasticity, P)
"""
function max_state_change(rs::RateStateFriction, re::AbstractElasticity, P)
    # following: Lapusta, N., Rice, J.R., Ben‐Zion, Y., Zheng, G., 2000. Elastodynamic analysis for slow tectonic loading with spontaneous rupture episodes on faults with rate‐ and state‐dependent friction. J. Geophys. Res. 105, 23765–23789. https://doi.org/10.1029/2000JB900250

    Peff = P * (1 - rs.λ)
    faultwidth = rs.D
    
    if isa(re, IncompressibleElasticity)
        G = re.G
        ν = 0.5
    else
        G = re.G
        ν = re.ν
    end

    k = 2 / π * (G / (1 - ν)) / faultwidth
    xi = 0.25 * (k * rs.L / (rs.a * Peff) - (rs.b - rs.a) / rs.a)^2 - k * rs.L / (rs.a * Peff)
    if xi > 0
        θmax = rs.a * Peff / (k * rs.L - (rs.b - rs.a) * Peff)
    elseif xi < 0
        θmax = 1 - (rs.b - rs.a) * Peff / (k * rs.L)
    else
        error("xi=0 in time step calculation")
    end
    θmax = clamp(θmax, 0.1, 0.2) # limit maximum according to Lapusta and Liu (2009); lower limit is not mentioned in the paper but it's in the code
    return θmax
end


function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1e6, dt_max=1e6, dt_min=1e-4)

    # Extract elastic stresses/pressure from solution vector
    τ    = zeros(ntime)
    τy   = zeros(ntime)  
    Ω_v  = zeros(ntime)
    t_v  = zeros(ntime)
    V_v  = zeros(ntime)
    τ_e  = others.τ0
    τ[1]  = τ_e[1]
    τy[1] = others.P[1] .* (c[1].μ₀ .+ c[1].b * 0)
    P_e   = others.P
    Ω   = 0.0;
    t    = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, Ω_old=Ω, P0 = P_e, P=(50e6,))       # other non-differentiable variables needed to evaluate the state functions

        # solve
        x = solve(c, x, vars, others, verbose = true, xnorm=xnorm)
        t += others.dt

        ε_ratestate = compute_strain_rate(c[1]; τ = x[1], Ω_old=others.Ω_old, P=others.P[1], dt=others.dt)
        τ_e = compute_stress_elastic(c, x, others)
        Vp = 2.0 * c[1].D * ε_ratestate                 # slip velocity
        
        Ω =  update_Ω(c[1]; ε=ε_ratestate, others...)   # update state variable

        # interface strength (for visualization)
        istrength = others.P[1] .* (c[1].μ₀ .+ c[1].b * Ω)

        τy[i]  = istrength
        τ[i]   = x[1]        # stress of the rate element
        Ω_v[i] = Ω
        t_v[i] = t
        V_v[i] = Vp

        # update time step (which depends on the state evolution)
        dt = compute_dt_ratestate(c[1], c[2], others.P[1], Vp, Ω, dt_min=1e-4, dt_max=1e6)     
    end

    return t_v, τ, Ω_v, V_v, τy
end

c, x, xnorm, vars, args, others = let

    a   = 0.011     # rate paramete
    b   = 0.015     # state paramter
    L   = 0.0047    # characteristic slip distance [m]
    μ₀  = 0.5       # reference friction
    V₀  = 4e-9      # reference slip rate [m/s]
    λ   = 0         # fluid pressure ratio
    C   = 0         # cohesion [Pa]
    D   = 500       # cell size

    hypoplasticity   = RateStateFriction(λ, μ₀, V₀, a, b, L, C, D)
    elastic          = IncompressibleElasticity(20e9)

    # Maxwell viscoelastic model
    # elastic --- viscous/rate_state
    c  = SeriesModel(hypoplasticity, elastic)
    
    platerate = 4e-9

    vars   = (; ε = platerate/D/2, θ = 1.0e-20)          # input variables (constant)
    args   = (; τ = 2.0e1, P = 50.0e6)                    # guess variables (we solve for these, differentiable)
    others = (; dt = 1.1e-2, τ0 = (2e7, ), P0 = (0.0, ), P=(50e6,)) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)
    char_τ  = 1e6
    char_ε  = vars.ε+vars.θ
    xnorm   = normalisation_x(c, char_τ, char_ε)

    c, x, xnorm, vars, args, others

end


t_v, τ, Ω_v, V_v, τy = stress_time(c, vars, x, xnorm, others; ntime = 1000, dt=1e6)


# Plotting
SecYear = 3600 * 24 * 365.25
fig = Figure(fontsize = 30, size = (1200, 600) .* 1)
ax1 = Axis(fig[1, 1], title = "Rate & state hypoplasticity", xlabel = L"$t$ [yr]", ylabel = L"$\tau_{II}$ [MPa]")
lines!(ax1, t_v/SecYear , τy/1e6, color=:red, label=L"\tau_{interface}")
lines!(ax1, t_v/SecYear , τ/1e6, color=:blue, label=L"\tau_{II}")
scatter!(ax1, t_v/SecYear , τ/1e6, color=:blue)
axislegend(ax1, position = :rb, labelsize=18)

ax2  = Axis(fig[1, 2],  xlabel = L"$t$ [yr]", ylabel = L"$\Omega$")
ax2r = Axis(fig[1, 2]; ylabel = L"$V_p$ [m/s]", yaxisposition = :right,
            ylabelcolor = :green, yticklabelcolor = :green, ytickcolor = :green, yscale = log10)
lines!(ax2,  t_v[2:end]/SecYear , Ω_v[2:end], color=:black, label = "omega")
lines!(ax2r, t_v[2:end]/SecYear , V_v[2:end], color=:green, label = "slip velocity")

display(fig)
