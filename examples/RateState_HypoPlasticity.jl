using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie
import Statistics: mean

include("../rheologies/RheologyDefinitions.jl")
include("../rheologies/RateState_HypoPlastic.jl")

function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1e-2)
    # Extract elastic stresses/pressure from solution vector
    τ1   = zeros(ntime)
    Ω_v  = zeros(ntime)
    t_v  = zeros(ntime)
    τ_e  = (25e3,)
    τ1[1] = τ_e[1]
    P_e  = (50e6,)
    Ω   = 0.0;
    t    = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, Ω_old=Ω, P0 = P_e, P=(50e6,))       # other non-differentiable variables needed to evaluate the state functions

        res = RheologyCalculator.compute_residual(c, x, vars, others)
        @show res

        x = solve(c, x, vars, others, verbose = true, xnorm=xnorm)
        τ1[i] = x[1]
        t += others.dt

        ε_ratestate = compute_strain_rate(c[1]; τ = x[1], Ω_old=others.Ω_old, P=others.P[1], dt=others.dt)
        τ_e = compute_stress_elastic(c, x, others)
        
        Ω =  update_Ω(c[1]; ε=ε_ratestate, others...) # update state variable
        Ω_v[i] = Ω
        t_v[i] = t
    end

    return t_v, τ1, Ω_v
end

c, x, xnorm, vars, args, others = let

    a   = 0.011 # rate paramete
    b   = 0.015 # state paramter
    L   = 0.0047 # characteristic slip distance [m]
    μ₀  = 0.5 # reference friction
    V₀  = 4e-9 # reference slip rate [m/s]
    λ   = 0 # fluid pressure ratio
    C   = 0 # cohesion [Pa]
    D   = 500   # cell size

    hypoplasticity   = RateStateFriction(λ, μ₀, V₀, a, b, L, C, D)
    elastic          = IncompressibleElasticity(20e9)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c  = SeriesModel(hypoplasticity, elastic)
    #c  = SeriesModel(elastic)
    #c  = SeriesModel(hypoplasticity)
    
    platerate = 4e-9

    vars   = (; ε = platerate/D/2, θ = 1.0e-20)          # input variables (constant)
    args   = (; τ = 2.0e1, P = 1.0e6)                    # guess variables (we solve for these, differentiable)
    others = (; dt = 1.1e-2, τ0 = (0e0, ), P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)
    char_τ  = 1e6
    char_ε  = vars.ε+vars.θ
    xnorm   = normalisation_x(c, char_τ, char_ε)

    c, x, xnorm, vars, args, others

end


t_v, τ, Ω_v = stress_time(c, vars, x, xnorm, others; ntime = 1000, dt=1e6)


SecYear = 3600 * 24 * 365.25
fig = Figure(fontsize = 30, size = (800, 600) .* 1)
ax1 = Axis(fig[1, 1], title = "Rate & state hypoplasticity", xlabel = L"$t$ [s]", ylabel = L"$\tau$ [MPa]")
lines!(ax1, t_v , τ/1e6, color=:black, label = "stress")
axislegend(ax1, position = :rb, labelsize=18)

ax2 = Axis(fig[1, 2],  xlabel = L"$t$ [s]", ylabel = L"$\Omega$")
lines!(ax2, t_v , Ω_v, color=:black, label = "omega")
axislegend(ax2, position = :rb, labelsize=18)

display(fig)


    