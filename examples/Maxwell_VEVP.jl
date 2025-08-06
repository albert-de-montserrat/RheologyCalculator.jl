using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import RheologyCalculator as RC

using GLMakie

analytical_solution(ϵ, t, G, η) = 2 * ϵ * η * (1 - exp(-G * t / η))

function stress_time(c, vars, x; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solutio vector
    τ1   = zeros(ntime)
    λ    = zeros(ntime)
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
        τ1[i]  = x[1]
        t += others.dt
        τ_an[i] = analytical_solution(vars.ε, t, c.leafs[2].G, c.leafs[1].η)
        τ_e     = x[1] # compute_stress_elastic(c, x, others)
    
        t_v[i] = t
    end
    return t_v, τ1, τ_an
end

c, x, vars, args, others = let
    viscous     = LinearViscosity(1e22)
    viscous_reg = LinearViscosity(1e20)
    elastic     = IncompressibleElasticity(10e9)
    plastic     = DruckerPrager(10e6, 30, 0)

    # Maxwell visco-elasto-(visco-plastic) model
    p  = ParallelModel(viscous_reg, plastic)
    c  = SeriesModel(viscous, elastic, p)
    # c  = SeriesModel(viscous, elastic)

    # input variables (constant)
    vars   = (; ε = 1.0e-14, θ = 1.0e-20)
    # guess variables (we solve for these, differentiable)
    args   = (; τ = 2.0e3, P = 1.0e6, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e8, τ0 = (0e0, ), P0 = (0.0, ))

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

using StaticArrays

# let
    x = SA[args.τ, vars.ε, 0.0, 0.0] # initial guess for the solver
    t_v, τ, τ_an = stress_time(c, vars, x; ntime = 1_000, dt = 1e8)

    SecYear = 3600 * 24 * 365.25
    fig = Figure(fontsize = 30, size = (800, 600) .* 2)
    ax  = Axis(fig[1, 1], title = "Visco-elasto-plastic model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")

    lines!(ax, t_v / SecYear / 1.0e3, τ_an / 1.0e6, color=:black, label = "viscoelastic analytical")
    scatter!(ax, t_v / SecYear / 1.0e3, τ / 1.0e6,  color=:red, label = "numerical")


    axislegend(ax, position = :rb)
    ax.xlabel = L"t [kyr]"
    ax.ylabel = L"\tau [MPa]"
    display(fig)
# # end

# using ForwardDiff

r = RC.compute_residual(c, x, vars, others)
# J = ForwardDiff.jacobian(y -> RC.compute_residual(c, y, vars, others), x)

# # RC.global_series_functions(c)
# # x
# eqs      = RC.generate_equations(c)
# args_all = RC.generate_args_template(eqs, x, others)

# args_all[1]
# args_all[2]
# args_all[3]
# args_all[4]
# x
# # evaluates the self-components of the residual
# residual1 = RC.evaluate_state_functions(eqs, args_all, others)
# residual2 = RC.add_children(residual1, x, eqs)
# residual3 = RC.subtract_parent(residual2, x, eqs, vars)


# add_child(x, eqs, eqs[i].child)

# for eq in eqs
#     println("
#     Equation: $(eq.fn)
#     Children: $(eq.child)
#     Parents: $(eq.parent)
#     "
#     )
# end

# RC.add_children(residual1, x, eqs)


# # @inline get_own_functions(c::ParallelModel) = get_own_functions(c, parallel_state_functions, global_parallel_state_functions, local_parallel_state_functions)

# # fns_own_all    = RC.parallel_state_functions(p.leafs)
# # fns_own_global = RC.global_parallel_state_functions(fns_own_all) |> RC.superflatten |> RC.flatten_repeated_functions
# # fns_own_local  = RC.local_parallel_state_functions(fns_own_all)
 

# #  v = Base.@ntuple $NC i -> begin
# i = 1
# eq_ind = child[i]
# RC.add_child(x, eqs[eq_ind], eq_ind)
# @edit RC.add_child(x, eqs[eq_ind], eq_ind)

# eqs[eq_ind] isa RC.CompositeEquation{A, B, typeof(RheologyCalculator.compute_lambda)} where {A, B}

# # end
# @show v, child
# sum(v)