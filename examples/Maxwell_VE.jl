using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("RheologyDefinitions.jl")

using GLMakie

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solutio vector
    τ1   = zeros(ntime)
    τ_an = zeros(ntime)
    # τ2 = zeros(ntime)
    # P1 = zeros(ntime)
    # P2 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (0.0,)
    P_e = (0.0,)
    t = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

        x = solve(c, x, vars, others, verbose = false)
        τ1[i] = x[1]
        t += others.dt
        τ_an[i] = analytical_solution(vars.ε, t, c.leafs[2].G, c.leafs[1].η)
        τ_e = compute_stress_elastic(c, x, others)
    
        t_v[i] = t
    end

    return t_v, τ1, τ_an
end

c, x, vars, args, others = let

    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c  = SeriesModel(viscous, elastic)

    vars = (; ε = 1.0e-14, θ = 1.0e-20)         # input variables (constant)
    args = (; τ = 2.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0e0, ), P0 = (0.0, ))       # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

t_v, τ, τ_an = stress_time(c, vars, x; ntime = 10_000, dt = 1e9)

SecYear = 3600 * 24 * 365.25
fig = Figure(fontsize = 30, size = (800, 600) .* 2)
ax  = Axis(fig[1, 1], title = "Maxwell VE model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")
ax2 = Axis(fig[2, 1], title = "Maxwell VE model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")

lines!(ax, t_v / SecYear / 1.0e3, τ_an / 1.0e6, color=:black, label = "analytical")
scatter!(ax, t_v / SecYear / 1.0e3, τ / 1.0e6,  color=:red, label = "numerical")

lines!(ax2, t_v / SecYear / 1.0e3, log10.(abs.(τ_an.-τ) ./ τ_an), color=:black)

axislegend(ax, position = :rb)
#title!(ax,"Burgers model")
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [MPa]"

ax2.xlabel = L"t [kyr]"
ax2.ylabel = L"\log_{10}\text{relative error}"
display(fig)