using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solutio vector
    τ1 = zeros(ntime)
    # τ2 = zeros(ntime)
    # P1 = zeros(ntime)
    # P2 = zeros(ntime)
    t_v = zeros(ntime)
    τ_e = (0.0, 0.0)
    # P_e = (0.0, 0.0)
    t = 0.0
    for i in 2:ntime
        #global τ_e, P_e, x, vars, others, t
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions

        x = solve(c, x, vars, others, verbose = false)
        τ1[i] = x[1]

        τ_e = compute_stress_elastic(c, x, others)
        P_e = compute_pressure_elastic(c, x, others)

        # τ1[i] = τ_e[1]
        # P1[i] = P_e[1]
        # if length(τ_e) > 1
        #     τ2[i] = τ_e[2]
        #     P2[i] = P_e[2]
        # end
        t += others.dt
        t_v[i] = t
    end

    return t_v, τ1
end

c, x, vars, args, others = let

    viscous1 = LinearViscosity(1e21)
    elastic0 = Elasticity(1.0e10, 4.0e10)
    elastic1 = Elasticity(1.0e10, 3.0e10)

    # Standard linear model (SML)
    #      parallel
    #          |
    # elastic --- viscous
    #    |
    # viscous

    c0 = SeriesModel(viscous1, elastic1)
    p  = ParallelModel(c0, elastic0)
    c  = SeriesModel(p)

    vars = (; ε = 1.0e-15, θ = 1.0e-20)         # input variables (constant)
    args = (; τ = 2.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (1.0, 2.0), P0 = (0.0, 0.1))       # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

t_v, τ1 = stress_time(c, vars, x; ntime = 2_000, dt = 1.0e8)


SecYear = 3600 * 24 * 365.25
fig = Figure(fontsize = 30, size = (800, 600) .* 2)
ax = Axis(fig[1, 1], title = "Burgers model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]", xticks = (0:2:20))

scatter!(ax, t_v / SecYear / 1.0e3, τ1 / 1.0e6, label = "τ1")
#scatter!(ax,t_v/SecYear/1e3,τ2/1e6, label="τ2")
# lines!(ax, t_anal / SecYear / 1.0e3, τ1_anal / 1.0e6, label = "τ1 analytical")

axislegend(ax, position = :rb)
#title!(ax,"Burgers model")
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [Pa]"
display(fig)