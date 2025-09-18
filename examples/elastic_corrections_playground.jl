# tests the elastic effective strainrate approach by computing the full tensor (in 2D, assuming that xx = -zz)
# and comparing it to invariant formulation
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie
import Statistics: mean
using LinearAlgebra

include("../rheologies/RheologyDefinitions.jl")

# Define rheology

c, x, τ0, ε, args, others = let

    viscous = LinearViscosity(1e21)
    elastic = IncompressibleElasticity(10e9)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c      = SeriesModel(viscous, elastic)

    vars   = (; ε = 1.0e-14, θ = 0.0e-20)                # input variables (constant)
    args   = (; τ = 2.0e3, P = 1.0e6)                    # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0e0, ), P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions


    x = initial_guess_x(c, vars, args, others)
    τ0 = ([1.0, 1e3],)
    ε  = [1e-14, 1e-14]
    c, x, τ0, ε, args, others
end

#=
c, x, τ0, ε, args, others = let

    viscous = LinearViscosity(1e20)
    elastic = IncompressibleElasticity(5e10)

    # Maxwell viscoelastic model
    # elastic --- viscous

    c      = SeriesModel(ParallelModel(viscous, elastic))   # 1 elastic kelvin element
    #c      = SeriesModel(viscous,elastic)

    vars   = (; ε = 1.0e-14, θ = 0.0e-20)                # input variables (constant)
    args   = (; τ = 2.0e3, P = 1.0e6)                    # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0e0, ), P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions


    x = initial_guess_x(c, vars, args, others)
    τ0 = ([1.0, 1e3],)
    ε  = [1e-14, 1e-14]
    c, x, τ0, ε, args, others
end
=#

second_invariant(a,b) = sqrt.(a.^2 .+ b.^2)

function stress_time_full_tensor(c, x, τ0, ε; ntime = 200, dt = 1.0e8)
    # This computes the full deviatoric stress tensor components (in 2D, assuming τzz = -τxx)
    # 
    τxx  = zeros(ntime) # note, we assume  τzz = -τxx
    τxz  = zeros(ntime)
    t_v  = zeros(ntime)

    τxx_e  = (τ0[1][1],)
    τxz_e  = (τ0[1][2],) 

    τxx[1] = τ0[1][1]
    τxz[1] = τ0[1][2]
    P_e  = (0.0,)
    t    = 0.0
    for i in 2:ntime
        others  = (; dt = dt, τ0 = τxx_e[1], P0 = P_e)              # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = ε[1], θ = 0.0) 
        x       = solve(c, x, vars, others, verbose = false, elastic_correction=false)
        τxx_e   = compute_stress_elastic(c, x, others)              # elastic stress components
        τxx[i]  = x[1]                                              # total stress    

        others  = (; dt = dt, τ0 = τxz_e[1], P0 = P_e)              # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = ε[2], θ = 0.0) 
        x       = solve(c, x, vars, others, verbose = false, elastic_correction=false)
        τxz_e   = compute_stress_elastic(c, x, others)              # elastic stress components
        τxz[i]  = x[1]                                              # total stress    


        t += others.dt

        t_v[i] = t
    end

    return t_v, τxx, τxz
end

function stress_time_invariant_manual(c, x0, τ0, ε; ntime = 200, dt = 1.0e8)
    # This computes the stress with invariants
    # 
    τxx  = zeros(ntime) # note, we assume  τzz = -τxx
    τxz  = zeros(ntime)
    τII  = zeros(ntime)
    t_v  = zeros(ntime)

    G = c[2].G

    τxx_e  = (τ0[1][1],)
    τxz_e  = (τ0[1][2],) 

    τxx[1] = τ0[1][1]
    τxz[1] = τ0[1][2]
    P_e  = (0.0,)
    t    = 0.0
    for i in 2:ntime
        εeff_xx = ε[1] + τxx_e[1]/(2*G*dt)
        εeff_xz = ε[2] + τxz_e[1]/(2*G*dt)
        εeff_II = second_invariant(εeff_xx, εeff_xz)
        nij     = [εeff_xx, εeff_xz]./εeff_II                  # normalized deviatoric direction tensor

        others  = (; dt = dt, τ0 = 0.0, P0 = P_e)              # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = εeff_II, θ = 0.0) 
        x       = solve(c, x0, vars, others, verbose = false, elastic_correction=false)

        τII_e   = compute_stress_elastic(c, x, others)          # elastic stress components
        
        τxx_e   = nij[1]*τII_e[1]                               # get deviatoric stresses from invariant
        τxz_e   = nij[2]*τII_e[1]
        
        τII[i]  = x[1]                                           # total stress    

        t += others.dt

        t_v[i] = t
    end

    return t_v, τII
end


function analytics_kelvin(t_v, c, others, τ0, ε, Δt)
    # 1 elastic kelvin element:
    #|--⟦▪̲̅▫̲̅▫̲̅▫̲̅¹--|
    #|--/\/\/¹--|    

    η = c.branches[1][1].η
    G = c.branches[1][2].G
    
    τvis_xx = 2η*ε[1]
    τvis_xz = 2η*ε[2]
    τel_xx  = 2G*t_v.*ε[1] .+ τ0[1][1]
    τel_xz  = 2G*t_v.*ε[2] .+ τ0[1][2]

    τxx   = τel_xx .+ τvis_xx
    τxz   = τel_xz .+ τvis_xz

    return second_invariant(τxx, τxz)
end

function analytics_maxwell(t_v, c, others, τ0, ε, Δt)
    # 1 elastic maxwell element:
    # --⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--

    η = c[1].η
    G = c[2].G

    τxx = zero(t_v)
    τxz = zero(t_v)
    τxx[1] = τ0[1][1]
    τxz[1] = τ0[1][2]
    for i = 2:length(t_v)
        τxx_o = τxx[i-1]
        τxz_o = τxz[i-1]
        
        τxx[i]  = inv(inv(2G*Δt) + inv(2η))*(ε[1] + τxx_o/(2G*Δt))
        τxz[i]  = inv(inv(2G*Δt) + inv(2η))*(ε[2] + τxz_o/(2G*Δt))
        
    end
  
    return second_invariant(τxx, τxz)
end

SecYear = 3600*24*365.25
dt = SecYear*10
t_v, τxx, τxz = stress_time_full_tensor(c, x, τ0, ε; ntime = 20, dt = dt) 
t_v, τII_invariants = stress_time_invariant_manual(c, x, τ0, ε; ntime = 20, dt = dt)
τII_tot = second_invariant(τxx, τxz)
#τII_ana = analytics_kelvin(t_v, c, others, τ0, ε, dt)
τII_ana = analytics_maxwell(t_v, c, others, τ0, ε, dt)

error  = norm(τII_tot[2:end] .- τII_ana[2:end])
errorI = norm(τII_invariants[2:end] .- τII_ana[2:end])

@info "FullTensor - Analytical:" error/mean(τII_tot)
@info "Invariant  - Analytical:" errorI/mean(τII_tot)

fig, ax, li = scatter(t_v ./ SecYear, τII_tot ./ 1e6, color=:blue, label="numerics, full tensor", markersize=10)
lines!(ax, t_v[1:end]/SecYear, τII_ana[1:end]/1e6, color=:red, label="analytics", linewidth=3)
lines!(ax, t_v[1:end]/SecYear, τII_invariants[1:end]/1e6, color=:green, label="invariants", linewidth=3)

ax.xlabel= "Time (years)"
ax.ylabel= "Second invariant (MPa)"
display(fig)