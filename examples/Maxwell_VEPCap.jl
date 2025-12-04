# This implements the mode1/mode2 plasticity model proposed in:
# Popov, A. A., Berlie, N., and Kaus, B. J. P.: A dilatant visco-elasto-viscoplasticity model with globally continuous tensile cap: 
#   stable two-field mixed formulation, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-2469, 2025.
using Test, LinearAlgebra
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using GLMakie
using StaticArrays

include("../rheologies/RheologyDefinitions.jl")
include("../rheologies/DruckerPragerCap.jl")

function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
    # Extract elastic stresses/pressure from solution vector
    τ1      = zeros(ntime)
    λ       = zeros(ntime)
    P1      = zeros(ntime)
    F1      = zeros(ntime)
    
    mode2   = zeros(ntime)
    t_v     = zeros(ntime)
    τ_e     = (0.0,)
    P_e     = (0.0e6,)
    P1[1]   = P_e[1]
    τ1[1]   = τ_e[1]
    x       = SA[τ1[1],0, P1[1]]
    t       = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
        
        x = RheologyCalculator.solve(c, x, vars, others, verbose = true, xnorm=xnorm)
        
        t += others.dt
        
        τ_e = compute_stress_elastic(c, x, others)
        P_e = compute_pressure_elastic(c, x, others)
        τ1[i] = τ_e[1]
        P1[i] = P_e[1]
        F1[i] = compute_F(c.leafs[end], τ1[i], P1[i])   
#        mode2[i] = ismode2_yield(c.leafs[end], τ_e[1], P_e[1])

        t_v[i] = t
    end

    return t_v, τ1, P1, F1, mode2
end

c, x, xnorm, vars, args, others = let

    viscous = LinearViscosity(1e23)
    #elastic = IncompressibleElasticity(10e9)
    elastic = Elasticity(1e10, 2e11)
    plastic = DruckerPrager(1e6, 30, 10)
    plastic = DruckerPragerCap(; C=1e6, ϕ=30.0, ψ=10.0, η_vp=0.0, Pt=-5e5) 

    # Maxwell viscoelastic model
    # elastic --- viscous

    c  = SeriesModel(viscous, elastic, plastic)
    #c  = SeriesModel(elastic, plastic)
    #c  = SeriesModel(elastic)
    
    # input variables (constant)
    vars = (; ε = 0*7.0e-14, θ = 7.0e-15)
    # guess variables (we solve for these, differentiable)
    args = (; τ = 0.0e3, P = 0.3e6, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e5, τ0 = (0e0, ), P0 = (0.3e6, ))

    x       = initial_guess_x(c, vars, args, others)
    char_τ  = plastic.C
    char_ε  = vars.ε+vars.θ
    xnorm   = normalisation_x(c, char_τ, char_ε)

    c, x, xnorm, vars, args, others
end



# Plot yield stress - this is reproducing Fig. 2 of the paper
τII = 0:0.01e6:2e6
P   = -1e6:0.01e6:2e6

F = zeros(length(P), length(τII))
Q = zeros(length(P), length(τII))
ismode2 = zeros(length(P), length(τII))
for i in CartesianIndices(F)
    F[i] = compute_F(c.leafs[end], τII[i[2]], P[i[1]])
    Q[i] = compute_Q(c.leafs[end], τII[i[2]], P[i[1]])
#    ismode2[i] = ismode2_yield(c.leafs[end], τII[i[2]], P[i[1]])
end
ismode2[F.>0] .= NaN
Q1 = copy(Q)
Q1[F.<1e-8] .= NaN

function figure()
    # Reproduce Fig. 5 of Popov et al. (2025)
    fig = Figure(fontsize = 20, size = (800, 800) )
    ax1 = Axis(fig[1,1], title="Volumetric extension (1)",  xlabel=L"$t$ [yr]",  ylabel=L"$P$, $\tau_{II}$ [MPa]", xlabelsize=20, ylabelsize=20)
    ax2 = Axis(fig[1,2], title="Volumetric compaction (2)", xlabel=L"$t$ [yr]",  ylabel=L"$P$, $\tau_{II}$ [MPa]", xlabelsize=20, ylabelsize=20)
    ax3 = Axis(fig[2,1], title="Vol. + Dev. shear (3)",     xlabel=L"$t$ [yr]",  ylabel=L"$P$, $\tau_{II}$ [MPa]", xlabelsize=20, ylabelsize=20)
    ax4 = Axis(fig[2,2], title="Stress paths",               xlabel=L"$P$ [MPa]", ylabel=L"$\tau_{II}$ [MPa]",      xlabelsize=20, ylabelsize=20)

    SecYear = 3600 * 24 * 365.25
    t_v1, τ1, P1, F1, mode2_1 = stress_time(c, (; ε = 0*7.0e-14, θ =   7.0e-15), x, xnorm, others; ntime = 11, dt = SecYear*2)
    println("-------")
    t_v2, τ2, P2, F2, mode2_2 = stress_time(c, (; ε =   7.0e-14, θ = 0*7.0e-15), x, xnorm, others; ntime = 80, dt = 1e7)
    println("-------")
    t_v3, τ3, P3, F3, mode2_3 = stress_time(c, (; ε =   7.0e-14, θ =   7.0e-15), x, xnorm, others; ntime = 30, dt = 2e7)
    println("-------")

    lines!(ax1, t_v1 / SecYear , P1 / 1.0e6, color=:red, label = L"$P$")
    lines!(ax1, t_v1 / SecYear , τ1 / 1.0e6,  color=:blue, label =  L"$\tau_{II}$")
    axislegend(ax1, position=:lb)

    lines!(ax2, t_v2 / SecYear , P2 / 1.0e6, color=:red, label = L"$P$")
    lines!(ax2, t_v2 / SecYear , τ2 / 1.0e6,  color=:blue, label =  L"$\tau_{II}$")

    lines!(ax3, t_v3 / SecYear , P3 / 1.0e6, color=:red, label = L"$P$")
    lines!(ax3, t_v3 / SecYear , τ3 / 1.0e6,  color=:blue, label =  L"$\tau_{II}$")

    GLMakie.contour!(ax4, P/1e6, τII/1e6, F, levels = [0.01], color = :black)
    GLMakie.contour!(ax4, P/1e6, τII/1e6, Q, levels = [0.001], color = :black, linestyle=:dash)
    GLMakie.scatter!(ax4, P1/1e6, τ1/1e6, color = :green, label=L"1")
    GLMakie.scatter!(ax4, P2/1e6, τ2/1e6, color = :red, label=L"2")
    GLMakie.scatter!(ax4, P3/1e6, τ3/1e6, color = :blue, label=L"3")

    display(fig)
end

with_theme(figure, theme_latexfonts())

