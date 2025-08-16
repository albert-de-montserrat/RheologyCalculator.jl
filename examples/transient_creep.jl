using ForwardDiff, RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("RheologyDefinitions.jl")

using GLMakie

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solutio vector
    τ1 = zeros(ntime)
    τ2 = zeros(ntime)
    P1 = zeros(ntime)
    P2 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (0.0, 0.0)
    P_e = (0.0, 0.0)
    t = 0.0
    for i in 2:ntime
        #global τ_e, P_e, x, vars, others, t
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

        x = solve(c, x, vars, others)
        τ_e = compute_stress_elastic(c, x, others)
        P_e = compute_pressure_elastic(c, x, others)

        @inbounds τ1[i] = τ_e[1]
        @inbounds P1[i] = P_e[1]
        # if length(τ_e) > 1
            @inbounds τ2[i] = τ_e[2]
            @inbounds P2[i] = P_e[2]
        # end
        t += others.dt
        t_v[i] = t
    end
    return t_v, τ1, τ2, P1, P2, x
end


viscous1 = LinearViscosity(5.0e19)
viscous2 = LinearViscosity(1.0e20)
elastic1 = IncompressibleElasticity(1.0e10)
elastic2 = IncompressibleElasticity(1.0e10)

c, x, vars, args, others = let
    #         parallel        -      parallel
    #            |                      |
    #   elastic --- viscous    elastic --- viscous

    p1 = ParallelModel(viscous1, elastic1)
    p2 = ParallelModel(viscous2, elastic2)
    c = SeriesModel(p1,p2)

    vars = (; ε = 1.0e-15, θ = 1.0e-20)                             # input variables (constant)
    args = (; τ = 2.0e3, P = 1.0e6)                                 # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0.0, 0.0), P0 = (0.0, 0.1))       # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

# Burgers model, numerics
t_v, τ1, τ2, P1, P2, x1 = stress_time(c, vars, x; ntime = 25, dt = 1.0e9);
t_anal, τ1_anal, τ2_anal = simulate_series_Burgers_model(G1, η1, G2, η2, vars.ε, t_v[end], (t_v[2] - t_v[1]) / 10);

SecYear = 3600 * 24 * 365.25
fig = Figure(fontsize = 30, size = (800, 600))
ax = Axis(fig[1, 1], title = "Burgers model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")

lines!(ax, t_anal / SecYear / 1.0e3, τ1_anal / 1.0e6, label = "analytical", linewidth = 5, color = :black)
scatter!(ax, t_v / SecYear / 1.0e3, τ1 / 1.0e6, label = "numerical", color = :red, markersize = 15)
#scatter!(ax,t_v/SecYear/1e3,τ2/1e6, label="τ2")

axislegend(ax, position = :rb)
#title!(ax,"Burgers model")
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [MPa]"
save("Burgers_model.png", fig)
display(fig)