# using GLMakie
using ForwardDiff, StaticArrays, LinearAlgebra

import Base.IteratorsMD.flatten

using RheologyCalculator

# include("rheology_types.jl")
# include("state_functions.jl")
# include("composite.jl")
# include("kwargs.jl")
# include("recursion.jl")
# include("equations.jl")
# include("others.jl")
# include("post_calculations.jl")
# include("initial_guess.jl")
# include("local_solver.jl")
# include("../src/print_rheology.jl")

viscous1 = LinearViscosity(5.0e19)
viscous2 = LinearViscosity(1.0e20)
viscousbulk = BulkViscosity(1.0e18)
powerlaw = PowerLawViscosity(5.0e19, 3)
drucker = DruckerPrager(1.0e6, 10.0, 0.0)
elastic  = Elasticity(1.0e10, 1.0e12)
elastic1 = Elasticity(1.0e11, 1.0e13)

elasticbulk = BulkElasticity(1.0e10)
elasticinc = IncompressibleElasticity(1.0e10)

LTP = LTPViscosity(6.2e-13, 76, 1.8e9, 3.4e9)
diffusion = DiffusionCreep(1, 1, 1, 1.5e-3, 1, 1, 1)
dislocation = DislocationCreep(3.5, 1, 1.1e-16, 1, 1, 1)

c, x, vars, args, others = let
    # elastic - viscous -- parallel
    #                         |
    #                viscous --- viscous
    s1 = SeriesModel(viscous1, viscous2)
    p = ParallelModel(viscous1, viscous2)
    c = SeriesModel(elastic, viscous1, p)
    vars = (; ε = 1.0e-15, θ = 1.0e-20)      # input variables (constant)
    args = (; τ = 1.0e3, P = 1.0e6) # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)[1], # local  guess(es)
    ]

    c, x, vars, args, others
end

# ε = exp10.(LinRange(log10(1.0e-15), log10(1.0e-8), 50))
# τ = similar(ε)
# x0 = copy(x)
# # args   = (; τ = 1e10)   # guess variables (we solve for these, differentiable)
# # others = (; dt = 1e10) # other non-differentiable variables needed to evaluate the state functions
# i=1
# # vars = (; ε = ε[i]) # input variables (constant)
# vars = (; ε = ε[i], θ = 1.0e-20) # input variables (constant)
# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    # elastic - viscous
    c = SeriesModel(elastic, viscous1)
    vars = (; ε = 1.0e-15, θ = 1.0e-20) # input variables (constant)
    args = (; τ = 1.0e2, P = 1.0e6)     # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10)            # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        # values(vars)..., # local  guess(es)
    ]

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    # viscous -- parallel
    #               |
    #      viscous --- viscous
    #         |
    #      viscous
    s1 = SeriesModel(viscous1, viscous2)
    p = ParallelModel(s1, viscous2)
    c = SeriesModel(viscous1, p)
    vars = (; ε = 1.0e-15) # input variables (constant)
    args = (; τ = 1.0e2) # guess variables (we solve for these, differentiable)
    others = (;)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
        values(args)..., # local  guess(es)
    ]

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    # viscous -- parallel
    #               |
    #      viscous --- viscous
    p = ParallelModel(viscous1, viscous2)
    c = SeriesModel(viscous1, p)
    vars = (; ε = 1.0e-15) # input variables (constant)
    args = (; τ = 1.0e2) # guess variables (we solve for these, differentiable)
    others = (;)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
    ]

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    # viscous -- parallel
    #               |
    #      viscous --- viscous
    #         |
    #      viscous
    s1 = SeriesModel(viscous1, viscous2)
    p = ParallelModel(s1, viscous2)
    c = SeriesModel(viscous1, p)
    vars = (; ε = 1.0e-15) # input variables (constant)
    args = (; τ = 1.0e2) # guess variables (we solve for these, differentiable)
    others = (;)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
        values(args)..., # local  guess(es)
    ]

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let # CRASHES
    #           parallel
    #               |
    #      viscous --- viscous
    #         |
    #      viscous
    s1 = SeriesModel(viscous1, viscous2)
    c = ParallelModel(s1, viscous2)

    vars = (; ε = 1.0e-15) # input variables (constant)
    args = (; τ = 1.0e2) # guess variables (we solve for these, differentiable)
    others = (;)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(vars)..., # local  guess(es)
        values(args)..., # global guess(es), solving for these
    ]

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    #           parallel
    #               |
    #      viscous --- viscous
    #         |
    #      viscous
    s1 = SeriesModel(viscous1, viscous2)
    c = ParallelModel(viscous1, s1) |> SeriesModel

    vars = (; ε = 1.0e-15) # input variables (constant)
    args = (; τ = 1.0e2) # guess variables (we solve for these, differentiable)
    others = (;)       # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    # viscous -- parallel    --      parallel
    #               |                   |
    #      viscous --- viscous  viscous --- viscous
    #         |
    #      viscous
    s1 = SeriesModel(viscous1, viscous2)
    p = ParallelModel(s1, viscous2)
    p1 = ParallelModel(viscous1, viscous2)
    c = SeriesModel(viscous1, p, p1)
    vars = (; ε = 1.0e-15) # input variables (constant)
    args = (; τ = 1.0e2) # guess variables (we solve for these, differentiable)
    others = (;)       # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)
    #x = SA[
    #    values(args)..., # global guess(es), solving for these
    #    values(vars)..., # local  guess(es)
    #    values(args)..., # local  guess(es)
    #    values(vars)..., # local  guess(es)
    #]

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)


#=
# Arne's model 1
c, x, vars, args, others = let
    # viscous -- parallel
    #               |  
    #      viscous --- viscous  
    #         |  
    #      viscous
    s1     = SeriesModel(viscous1, viscous2, viscous1)
    p      = ParallelModel(s1, viscous2)
    c      = SeriesModel(elastic, p)
    vars   = (; ε = 1e-15, θ = 1e-20) # input variables (constant)
    args   = (; τ = 1e2,   P = 1e6) # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10)       # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
        values(vars)..., # local  guess(es)
        values(args)..., # local  guess(es)
    ]

    c, x, vars, args, others
end

# Arne's model 2
c, x, vars, args, others = let
    viscous1    = LinearViscosity(1e12)
    viscous2    = LinearViscosity(1e12)
    elastic     = Elasticity(100e9, 1e12)
    LTP         = LTPViscosity(6.2e-13, 76, 1.8e9, 3.4e9)
    diffusion   = DiffusionCreep(1, 1, 1, 1.5e-3, 1, 1, 1)
    dislocation = DislocationCreep(3.5, 1, 1.1e-16, 1, 1, 1)

    # viscous -- parallel
    #               |  
    #      viscous --- viscous  
    #         |  
    #      viscous
    # s1     = SeriesModel(diffusion, LTP, dislocation)
    # p1     = ParallelModel(s1, viscous2)
    # p2     = ParallelModel(elastic, viscous2)
    # c      = SeriesModel(p1, p2)

    s1     = SeriesModel(diffusion, dislocation)
    p1     = ParallelModel(s1, viscous2)
    p2     = ParallelModel(elastic, viscous2)
    c      = SeriesModel(p1, p2)
    vars   = (; ε = 1e-12 * 2) # input variables (constant)
    args   = (; τ = 2e9)   # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10) # other non-differentiable variables needed to evaluate the state functions

    #x = SA[
    #    values(args)..., # global guess(es), solving for these
    #    0 .* values(vars)..., # local  guess(es)
    #    0 .* values(vars)..., # local  guess(es)
    #    0 .* values(args)..., # local  guess(es)
    #]
    x = initial_guess_x(c, vars, args, others)


    c, x, vars, args, others
end


# Arne's model 3
c, x, vars, args, others = let
    viscous1    = LinearViscosity(1e12)
    viscous2    = LinearViscosity(1e12)
    elastic     = Elasticity(100e9, 1e12)
    LTP         = LTPViscosity(6.2e-13, 76, 1.8e9, 3.4e9)
    diffusion   = DiffusionCreep(1, 1, 1, 1.5e-3, 1, 1, 1)
    dislocation = DislocationCreep(3.5, 1, 1.1e-16, 1, 1, 1)
    
    vars   = (; ε = 1e-12) # input variables (constant)
    args   = (; τ = 2e9)       # guess variables (we solve for these, differentiable)
    others = (; dt = 1e-2)     # other non-differentiable variables needed to evaluate the state functions

    # c      = SeriesModel(LTP)
    # x = SA[
    #     values(args)..., # global guess(es), solving for these
    # ]

    # p      = ParallelModel(LTP, viscous1)
    # c      = SeriesModel(p)
    # x = SA[
    #     values(args)..., # global guess(es), solving for these
    #     values(vars)..., # local  guess(es)
    # ]
    
    # s1     = SeriesModel(diffusion, dislocation, LTP)
    # p      = ParallelModel(s1, viscous1)
    # c      = SeriesModel(p)
    # x = SA[
    #     values(args)..., # global guess(es), solving for these
    #     1 .* values(vars)..., # local  guess(es)
    #     1 .* values(args)..., # global guess(es), solving for these
    # ]

    # s1     = SeriesModel(diffusion, dislocation, LTP)
    # p      = ParallelModel(s1, viscous1)
    # c      = SeriesModel(p, viscous1)
    # x = SA[
    #     values(args)..., # global guess(es), solving for these
    #     1 .* values(vars)..., # local  guess(es)
    #     1 .* values(args)..., # global guess(es), solving for these
    # ]
    
    s1     = SeriesModel(diffusion, dislocation, LTP)
    p1     = ParallelModel(s1, viscous2)
    p2     = ParallelModel(elastic, viscous2)
    c      = SeriesModel(p1, elastic)

    vars   = (; ε = 1e-12 * 2, θ = 1e-20)  # input variables (constant)
    args   = (; τ = 1e9, P = 1e6) # guess variables (we solve for these, differentiable)
    #args   = (; τ = 1e9, ) # guess variables (we solve for these, differentiable)
    
    others = (; dt = 1e2)        # other non-differentiable variables needed to evaluate the state functions

    # solution vector
    x = SA[
        values(args)..., # global guess(es), solving for these
        1 .* values(vars)[1]..., # local  guess(es)
        1 .* values(args)[1]..., # local  guess(es)
    ]

    c, x, vars, args, others
end


c, x, vars, args, others = let
    # elastic - viscous -- bulkviscous -- bulkelastic 
    c      = SeriesModel(viscous1, elastic, viscousbulk, elasticbulk)
    vars   = (; ε = 1e-15, θ = 1e-20)       # input variables (constant)
    args   = (; τ = 1e3,   P = 1e6)         # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10)                  # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
    ]

    c, x, vars, args, others
end

c, x, vars, args, others = let
    #             parallel    
    #                |       
    #   viscousbulk --- elasticbulkelastic 
    c      = SeriesModel(ParallelModel(viscousbulk, elasticbulk))
    vars   = (; θ = 1e-20)       # input variables (constant)
    args   = (; P = 1e6)         # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10)                  # other non-differentiable variables needed to evaluate the state functions

    x = SA[
        values(args)..., # global guess(es), solving for these
    ]

    c, x, vars, args, others
end


c1, x, vars, args, others = let
    #      elastic - viscous -    parallel    
    #                                |       
    #                   viscousbulk --- elasticbulk

    p      = ParallelModel(viscousbulk, elasticbulk)
    c      = SeriesModel(viscous1, elastic, p)
    vars   = (; ε = 1e-15, θ = 1e-20)       # input variables (constant)
    args   = (; τ = 1e3,   P = 1e6)         # guess variables (we solve for these, differentiable)
    others = (; dt = 1e10)                  # other non-differentiable variables needed to evaluate the state functions

    #x = initial_guess_x(c, others, args, vars)
    x = SA[
        values(args)..., # global guess(es), solving for these
    ]
    #x = init

    c, x, vars, args, others
end
=#


c, x, vars, args, others = let
    # Burger's model
    #      elastic - viscous -    parallel
    #                                |
    #                   elastic --- viscous

    p = ParallelModel(viscous2, elastic)
    viscous3 = LinearViscosity(1.0e21)
    c = SeriesModel(viscous3, elastic1, p)
    #c = SeriesModel(viscous3, elastic1)

    #c = p
   # c = SeriesModel(p);
    #c = SeriesModel(elastic, p)
    
    vars = (; ε = 1.0e-15, θ = 1.0e-20)         # input variables (constant)
    args = (; τ = 1.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (1.0, 2.0), P0=(0.0, 0.1))       # other non-differentiable variables needed to evaluate the state functions
    
    
    x   = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

c, x, vars, args, others = let
    # Linear viscelastoplastic model
    #      elastic - viscous -    drucker

    viscous3 = LinearViscosity(1.0e21)
    drucker1 = DruckerPrager(1.0e6, 0.0, 0.0)
    c = SeriesModel(viscous3, elastic, drucker1)
    #c = SeriesModel(viscous1, elastic)
    
    
    vars = (; ε = 1.0e-15, θ = 1.0e-20)         # input variables (constant)
    args = (; τ = 1.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = 1.0, P0=0.0)       # other non-differentiable variables needed to evaluate the state functions
    
    
    x   = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end

# @b solve($(c, x, vars, others)...)

ε = exp10.(LinRange(log10(1.0e-15), log10(1.0e-8), 50))
τ = similar(ε)
x0 = copy(x)
# args   = (; τ = 1e10)   # guess variables (we solve for these, differentiable)
# others = (; dt = 1e10) # other non-differentiable variables needed to evaluate the state functions
i=1
# vars = (; ε = ε[i]) # input variables (constant)
vars = (; ε = ε[i], θ = 1.0e-20) # input variables (constant)
sol = solve(c, x, vars, others)

function main(c, x, vars, args, others)
    ε = exp10.(LinRange(log10(1.0e-15), log10(1.0e-8), 50))
    τ = similar(ε)
    x0 = copy(x)
    # args   = (; τ = 1e10)   # guess variables (we solve for these, differentiable)
    # others = (; dt = 1e10) # other non-differentiable variables needed to evaluate the state functions
    for i in eachindex(ε)
        # vars = (; ε = ε[i]) # input variables (constant)
        vars = (; ε = ε[i], θ = 1.0e-20) # input variables (constant)
        sol = solve(c, x, vars, others)
        x = x0
        τ[i] = sol[1]
    end

    f, ax, h = scatterlines(log10.(ε), log10.(τ))
    # f,ax,h = scatterlines(ε, τ)
    # ax.xlabel = L"\dot\varepsilon_{II}"
    # ax.ylabel = L"\tau_{II}"
    return f
end


#main(c, x, vars, args, others)
#eqs = generate_equations(c)

#x0   = initial_guess_x(c, vars, args, others)
r    = compute_residual(c, x, vars, others)
#xnew = solve(c, x0, vars, others)

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


# Burgers model, numerics
t_v, τ1, τ2, P1, P2, x1 = stress_time( vars, others, x; ntime=200, dt=1e9)


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

SecYear =  3600*24*365.25
fig     =  Figure(fontsize=30)
ax      =  Axis(fig[1, 1], title="Burgers model", xlabel="t [kyr]", ylabel=L"\tau [MPa]", xticks=(0:2:20))  

scatter!(ax,t_v/SecYear/1e3,τ1/1e6, label="τ1")
#scatter!(ax,t_v/SecYear/1e3,τ2/1e6, label="τ2")
#lines!(ax,t_anal/SecYear/1e3,τ1_anal/1e6, label="τ1 analytical")
#lines!(ax,t_anal/SecYear/1e3,τ2_anal/1e6, label="τ2 analytical")

axislegend(ax, position=:rb)
#title!(ax,"Burgers model")
#ax.xlabel = L"t [kyr]"
#ax.ylabel = L"\tau [Pa]"
display(fig)


c, x, vars, args, others = let
    # Linear viscelastoplastic model
    #      elastic - viscous -    ParallelModel
    #                                |          
    #                        drucker --- viscous 
    viscous3 = LinearViscosity(1.0e21)
    drucker1 = DruckerPrager(1.0e6, 30.0, 10.0)
    p = ParallelModel(drucker1, viscous1)
    c = SeriesModel(viscous3, elastic, p)
    
    vars = (; ε = 1.0e-15, θ = 1.0e-18)         # input variables (constant)
    args = (; τ = 1.0e3, P = 1.0e6)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = 1.0, P0=0.0)       # other non-differentiable variables needed to evaluate the state functions
    
    
    x   = initial_guess_x(c, vars, args, others)

    c, x, vars, args, others
end
J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)