
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie

function stress_time( vars, others, x; ntime=200, dt=1e8)
    # Extract elastic stresses/pressure from solutio vector
    τ1 = zeros(ntime)
    τ2 = zeros(ntime)
    P1 = zeros(ntime)
    P2 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (0.0, 0.0)
    P_e = (0.0, 0.0)
    t = 0.0
    for i=2:ntime
        #global τ_e, P_e, x, vars, others, t
        others = (; dt = dt, τ0 = τ_e, P0=P_e)       # other non-differentiable variables needed to evaluate the state functions
        
        x   = solve(c, x, vars, others, verbose=false)
        τ_e = compute_stress_elastic(c, x, others)
        P_e = compute_pressure_elastic(c, x, others)

        τ1[i] = τ_e[1]
        P1[i] = P_e[1]
        if length(τ_e) > 1
            τ2[i] = τ_e[2]
            P2[i] = P_e[2]
        end
        t += others.dt
        t_v[i] = t
    end

    return t_v, τ1, τ2, P1, P2, x
end


# Analytical solution for Burgers model
function simulate_series_Burgers_model(E1, η1, E2, η2, ε̇, t_max, dt)
    N       = Int(round(t_max / dt)) + 1
    t       = range(0, step=dt, length=N)
    σ       = zeros(N)          # Stress
    ε_KV    = zeros(N)          # Strain in Kelvin–Voigt element
    σ_KV    = zeros(N) 
    σ_spring = zeros(N) # Stress in the spring of the Kelvin-Voigt element 
    for i in 2:N
        # Previous values
        ε_KV_prev   = ε_KV[i-1]
        σ_prev      = σ[i-1]
        dεKVdt      = (σ_prev - E2 * ε_KV_prev) / η2    # Kelvin–Voigt strain rate
        ε_KV[i]     = ε_KV_prev + dt * dεKVdt           # Update ε_KV
        dσdt        = E1 * (ε̇ - σ_prev / η1 - dεKVdt)   # Stress rate from Maxwell element
        σ[i]        = σ_prev + dt * dσdt                # Update stress
        σ_KV[i]     = E2 * ε_KV[i] + η2 * dεKVdt        # Calculate σ_KV explicitly at this step
        # Stress in the spring of the Kelvin-Voigt element (elastic part)
        σ_spring[i] = E2 * ε_KV[i]
    end

    return t, σ, σ_spring
end

viscous1 = LinearViscosity(5.0e19)
viscous2 = LinearViscosity(1.0e20)
elastic  = Elasticity(1.0e10, 3.0e10)
elastic1 = Elasticity(1.0e10, 4.0e10)

c, x, vars, args, others = let
    # Burger's model
    #      elastic - viscous -    parallel
    #                                |
    #                   elastic --- viscous

    p        = ParallelModel(viscous2, elastic)
    viscous3 = LinearViscosity(1.0e21)
    c        = SeriesModel(viscous3, elastic1, p)
    
    vars = (; ε = 1.0e-15, θ = 1.0e-20)         # input variables (constant)
    args = (; τ = 2.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (1.0, 2.0), P0=(0.0, 0.1))       # other non-differentiable variables needed to evaluate the state functions
    
    x   = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

η1 =  2*c.leafs[1].η
G1 =  2*c.leafs[2].G
if length(c.branches) > 0
    η2 =  2*c.branches[1].leafs[1].η
    G2 =  2*c.branches[1].leafs[2].G
else
    η2 =  0.0
    G2 =  0.0
end

t_anal,τ1_anal,τ2_anal  = simulate_series_Burgers_model(G1, η1, G2, η2, vars.ε, t_v[end], (t_v[2]-t_v[1])/10)


# Burgers model, numerics
t_v, τ1, τ2, P1, P2, x1 = stress_time( vars, others, x; ntime=200, dt=1e9)


SecYear =  3600*24*365.25
fig     =  Figure(fontsize=30, size=(800, 600) .* 2)
ax      =  Axis(fig[1, 1], title="Burgers model", xlabel="t [kyr]", ylabel=L"\tau [MPa]", xticks=(0:2:20))  

scatter!(ax,t_v/SecYear/1e3,τ1/1e6, label="τ1", color=:red)
#scatter!(ax,t_v/SecYear/1e3,τ2/1e6, label="τ2")
lines!(ax,t_anal/SecYear/1e3,τ1_anal/1e6, label="τ1 analytical")

axislegend(ax, position=:rb)
#title!(ax,"Burgers model")
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [Pa]"
display(fig)

import RheologyCalculator as RC

# function solve(c::AbstractCompositeModel, x::SVector, vars, others; tol::Float64 = 1.0e-9, itermax = 1.0e4, verbose::Bool = false)

    it = 0
    er = Inf
    # local α
    # while er > tol
    #     it += 1
        r = RC.compute_residual(c, x, vars, others)
        J = ForwardDiff.jacobian(y -> RC.compute_residual(c, y, vars, others), x)
        Δx = J \ r
        α = bt_line_search(Δx, J, x, r, c, vars, others)
        x -= α .* Δx
        # check convergence
        er = mynorm(Δx, x)

#         it > itermax && break
#     end
#     if verbose
#         println("Iterations: $it, Error: $er, α = $α")
#     end
#     return x
# end
