# tests the elastic effective strainrate approach by computing the full tensor (in 2D using 3 components)
# and comparing it to the analytical solution and an invariant formulation
#
# This example:
# Maxwell & Kelvin viscoelastic model
# --⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--|--⟦▪̲̅▫̲̅▫̲̅▫̲̅²--|
#                     |--/\/\/²--| 

using RheologyCalculator, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

using GLMakie
import Statistics: mean
using LinearAlgebra

include("../rheologies/RheologyDefinitions.jl")

# Define rheology - combined Kelvin & Maxwell
c, x, τ0, ε, args, others = let

    viscous  = LinearViscosity(1e31)
    elastic  = IncompressibleElasticity(10e9)
    viscous2 = LinearViscosity(1e18)
    elastic2 = IncompressibleElasticity(5e8)
    
    # Maxwell & Kelvin viscoelastic model
    # --⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--|--⟦▪̲̅▫̲̅▫̲̅▫̲̅²--|
    #                     |--/\/\/²--|

    p      = ParallelModel(viscous2, elastic2)
    c      = SeriesModel(viscous, elastic, p)

    vars   = (; ε = 1.0e-14, θ = 0.0e-20)                # input variables (constant)
    args   = (; τ = 2.0e3, P = 1.0e6)                    # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0e0, ), P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions

    x = initial_guess_x(c, vars, args, others)
    τ0 = ([0.0, 0e3, 0.0],[1.0e5, -1e5, 9.0e5])  # initial stress state (2D: xx,zz,xz)
    ε  = [1e-14, -1e-14, 1e-15]         # applied dev. strainrates (xx,zz,xz)
    c, x, τ0, ε, args, others
end


second_invariant(xx,xz)    = sqrt.(xx.^2 .+ xz.^2)
second_invariant(xx,zz,xz) = sqrt.(0.5.*xx.^2 .+ 0.5.*zz.^2 .+ xz.^2)

function stress_time_full_tensor(c, x, τ0, ε; ntime = 200, dt = 1.0e8)
    # This computes the full deviatoric stress tensor components (in 2D)
    τxx   = zeros(ntime) 
    τzz   = zeros(ntime) 
    τxz   = zeros(ntime)
    τxx1  = zeros(ntime) 
    τzz1  = zeros(ntime) 
    τxz1  = zeros(ntime)
    τxx2  = zeros(ntime) 
    τzz2  = zeros(ntime)
    τxz2  = zeros(ntime)
    
    t_v  = zeros(ntime)

    # initial elastic stresses of 1th element
    τxx_e  = (τ0[1][1],τ0[2][1])
    τzz_e  = (τ0[1][2],τ0[2][2])
    τxz_e  = (τ0[1][3],τ0[2][3]) 

    τxx1[1] = τ0[1][1]
    τzz1[1] = τ0[1][2]
    τxz1[1] = τ0[1][3]
    τxx2[1] = τ0[2][1]
    τzz2[1] = τ0[2][2]
    τxz2[1] = τ0[2][3]
    
    P_e  = (0.0,0.0)
    t    = 0.0
    for i in 2:ntime

        # xx
        others  = (; dt = dt, τ0 = τxx_e, P0 = P_e)              # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = ε[1], θ = 0.0) 
        x       = solve(c, x, vars, others, verbose = true, elastic_correction=false, atol=1e-15, rtol = 1e-11, itermax=100)
        τxx_e   = compute_stress_elastic(c, x, others)              # elastic stress components
        τxx[i]  = x[1]                                          # total stress    
        τxx1[i] = τxx_e[1]                                          # total stress    
        τxx2[i] = τxx_e[2]                                          # total stress    
        
        # zz
        others  = (; dt = dt, τ0 = τzz_e, P0 = P_e)              # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = ε[2], θ = 0.0) 
        x       = solve(c, x, vars, others, verbose = true, elastic_correction=false, atol=1e-15, rtol = 1e-11, itermax=100)
        τzz_e   = compute_stress_elastic(c, x, others)              # elastic stress components
        τzz[i]  = x[1]                                          # total stress    
        τzz1[i] = τzz_e[1]                                          # total stress    
        τzz2[i] = τzz_e[2]                                          # total stress    
        
        # xz
        others  = (; dt = dt, τ0 = τxz_e, P0 = P_e)              # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = ε[3], θ = 0.0) 
        x       = solve(c, x, vars, others, verbose = false, elastic_correction=false, atol=1e-15, rtol = 1e-11, itermax=100)
        τxz_e   = compute_stress_elastic(c, x, others)              # elastic stress components
        τxz[i]  = x[1]                                              # total stress    
        τxz1[i] = τxz_e[1]                                          # total stress    
        τxz2[i] = τxz_e[2]                                          # total stress    
      
        t += others.dt

        t_v[i] = t
    end

    return t_v, τxx, τzz, τxz, τxx1, τzz1, τxz1, τxx2, τzz2, τxz2
end

function stress_time_invariant_manual(c, x0, τ0, ε; ntime = 200, dt = 1.0e8)
    # This computes the stress with invariants for a series case
    # 
    τxx    = zeros(ntime) # note, we assume  τzz = -τxx
    τzz    = zeros(ntime) # note, we assume  τzz = -τxx
    τxz    = zeros(ntime)
    τII    = zeros(ntime)
    τII_2  = zeros(ntime)
    
    t_v     = zeros(ntime)

    η1     = c[1].η
    G1     = c[2].G
    η2     = c.branches[1][1].η
    G2     = c.branches[1][2].G

    τxx_e  = (τ0[1][1], τ0[2][1])
    τzz_e  = (τ0[1][2], τ0[2][2])
    τxz_e  = (τ0[1][3], τ0[2][3]) 
    τII_e  = second_invariant(τxx_e, τzz_e, τxz_e)
    τxx[1] = τ0[1][1]
    τzz[1] = τ0[1][2]
    τxz[1] = τ0[1][3]
    τII[1] = second_invariant(τxx[1], τzz[1], τxz[1])

    P_e  = (0.0,)
    t    = 0.0
    for i in 2:ntime
        εeff_xx = ε[1] + τxx_e[1]/(2*G1*dt)  + τxx_e[2]/(2*G2*dt + 2*η2)
        εeff_zz = ε[2] + τzz_e[1]/(2*G1*dt)  + τzz_e[2]/(2*G2*dt + 2*η2)
        εeff_xz = ε[3] + τxz_e[1]/(2*G1*dt)  + τxz_e[2]/(2*G2*dt + 2*η2)
        εeff_II = second_invariant(εeff_xx, εeff_zz, εeff_xz)
        εII     = second_invariant(ε...)
        nij     = [εeff_xx, εeff_zz, εeff_xz]./εeff_II           # normalized deviatoric direction tensor

        others  = (; dt = dt, τ0 = 0.0, P0 = P_e)                # other non-differentiable variables needed to evaluate the state functions
        vars    = (; ε = εeff_II, θ = 0.0) 
        x       = solve(c, x0, vars, others, verbose = false, elastic_correction=false)
        τII[i]  = x[1]                                      # total stress from invariant formulation
        
        
        τII_e   = compute_stress_elastic(c, x, others) 
       # @show ΔτII_e    
        #τII_e = τII_e .+ ΔτII_e
        @show τII_e, x[1]

        #=
        τ0 = second_invariant(τxx_e, τzz_e, τxz_e)
        @show τ0
        #vars    = (; ε = εeff_II, θ = 0.0) 
        others  = (; dt = dt, τ0 = second_invariant(τxx_e, τzz_e, τxz_e), P0 = P_e)      
        τII_e   = compute_stress_elastic(c, x, others)           # elastic stress components
        # this is not correct
        εeff_xx = ε[1] + τxx_e[1]/(2*G1*dt)  
        εeff_zz = ε[2] + τzz_e[1]/(2*G1*dt)  
        εeff_xz = ε[3] + τxz_e[1]/(2*G1*dt) 
        εeff_II = second_invariant(εeff_xx, εeff_zz, εeff_xz)
        =#
        #εeff_kelvin = second_invariant(ε...) - τII[i]/(2*η1) - (τII[i] - τII[i-1])/(2*G1*dt)    # effective strainrate in kelvin body
        #τII_e_2 = εeff_kelvin*2*G2*dt + second_invariant(τxx_e[2], τzz_e[2], τxz_e[2])          # elastic stress components from kelvin body

        εeff_kelvin = 
        #τII_e_2v2 = x[2]*2*G2*dt 
        #τII_e = (τII[i], τII_e_2)

        #@show εeff_kelvin, εeff_II, τII_e_2v2, τII_e_2
        τxx_e   = nij[1].*τII_e                                  # get deviatoric stresses from invariant
        τzz_e   = nij[2].*τII_e                                  # get deviatoric stresses from invariant
        τxz_e   = nij[3].*τII_e 

                                                  # total stress
        τII_2[i] = second_invariant(τxx_e[2], τzz_e[2], τxz_e[2])   # total stress from elastic stresses

        t += others.dt

        t_v[i] = t
    end

    return t_v, τII, τII_2
end

function analytics_maxwell_kelvin(t_v, c, others, τ0, ε, Δt)
    # 1 elastic maxwell element:
    # --⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--|--⟦▪̲̅▫̲̅▫̲̅▫̲̅²--|
    #                     |--/\/\/²--|
    #
    # 
    # in x, the following is valid:
    #   εxx = εxx_vis + εxx_elastic + εxx_kelvin
    #   εxx = τxx/(2η1) + (τxx-τxx_e_old[1])/(2G1Δt)  + εxx_kelvin
    #   τxx = (2G2*Δt + 2η2)*εxx_kelvin + τxx_e_old[2]   => εxx_kelvin = (τxx - τxx_e_old[2])/(2G2*Δt + 2η2)
    #
    # So:
    #   εxx = τxx/(2η1) + (τxx-τxx_e_old[1])/(2G1Δt)  + (τxx - τxx_e_old[2])/(2G2*Δt + 2η2)
    η1  = c[1].η
    G1  = c[2].G
    η2  = c.branches[1][1].η
    G2  = c.branches[1][2].G

    τxx = zero(t_v)
    τzz = zero(t_v)
    τxz = zero(t_v)
    τxx[1] = τ0[1][1]
    τzz[1] = τ0[1][2]
    τxz[1] = τ0[1][3]
    τxx2_o = τ0[2][1]
    τzz2_o = τ0[2][2]
    τxz2_o = τ0[2][3]   
    for i = 2:length(t_v)
        τxx_o = τxx[i-1]    # of series element
        τzz_o = τzz[i-1]    # of series element
        τxz_o = τxz[i-1]    # of series element

        # kelvin body:
        τxx[i]  = inv(inv(2G1*Δt) + inv(2η1) + inv(2G2*Δt + 2η2))*(ε[1] + τxx_o/(2G1*Δt) + τxx2_o/(2G2*Δt + 2η2))
        τzz[i]  = inv(inv(2G1*Δt) + inv(2η1) + inv(2G2*Δt + 2η2))*(ε[2] + τzz_o/(2G1*Δt) + τzz2_o/(2G2*Δt + 2η2))
        τxz[i]  = inv(inv(2G1*Δt) + inv(2η1) + inv(2G2*Δt + 2η2))*(ε[3] + τxz_o/(2G1*Δt) + τxz2_o/(2G2*Δt + 2η2))
        

        εxx_kelvin = ε[1] - (τxx[i] - τxx_o)/(2G1*Δt) - τxx[i]/2η1
        εzz_kelvin = ε[2] - (τzz[i] - τzz_o)/(2G1*Δt) - τzz[i]/2η1
        εxz_kelvin = ε[3] - (τxz[i] - τxz_o)/(2G1*Δt) - τxz[i]/2η1

        # update τxx2_o, τxz2_o
        τxx2_o = τxx2_o + 2G2*Δt*εxx_kelvin 
        τzz2_o = τzz2_o + 2G2*Δt*εzz_kelvin 
        τxz2_o = τxz2_o + 2G2*Δt*εxz_kelvin 
        
    end

    return second_invariant(τxx,τzz,τxz)
end

SecYear = 3600*24*365.25
dt = SecYear*100

SecYear = 3600*24*365.25
dt = SecYear*10
t_v, τxx, τzz, τxz, τxx1, τzz1, τxz1, τxx2, τzz2, τxz2 = stress_time_full_tensor(c, x, τ0, ε; ntime = 25, dt = dt) 
t_v1, τII_invariants, τII_invariants_2 = stress_time_invariant_manual(c, x, τ0, ε; ntime = 25, dt = dt)
#t_v, τII_invariants, _ = stress_time_invariant_manual_parallel(c, x, τ0, ε; ntime = 20, dt = dt)
τII_tot   = second_invariant(τxx, τzz, τxz)
τII_1_tot = second_invariant(τxx1, τzz1, τxz1)
τII_2_tot = second_invariant(τxx2, τzz2, τxz2)

println("-----")
τII_ana = analytics_maxwell_kelvin(t_v, c, others, τ0, ε, dt)
#t_v, τII_invariants = stress_time_invariant_manual_parallel(c, x, τ0, ε; ntime = 5, dt = dt)



#τII_ana = analytics_kelvin(t_v, c, others, τ0, ε, dt)

error  = norm(τII_tot[2:end] .- τII_ana[2:end])
@info "FullTensor - Analytical:" error/mean(τII_tot)
errorI = norm(τII_invariants[2:end] .- τII_ana[2:end])
@info "Invariant  - Analytical:" errorI/mean(τII_tot)
#=
errorII= norm(τII_invariants[2:end] .- τII_tot[2:end])

@info "Invariant  - FullTensor:" errorII/mean(τII_tot)
=#

fig, ax, li = scatter(t_v[2:end]/SecYear, τII_tot[2:end]/1e6, color=:blue, label="numerics, full tensor", linewidth=3)
lines!(ax, t_v[2:end]/SecYear, τII_1_tot[2:end]/1e6, color=:magenta, label="analytics, elastic1", linewidth=3)
lines!(ax, t_v[2:end]/SecYear, τII_2_tot[2:end]/1e6, color=:orange, label="analytics, elastic2", linewidth=3)

lines!(ax, t_v[2:end]/SecYear, τII_ana[2:end]/1e6, color=:red, label="analytics", linewidth=3)
lines!(ax, t_v[2:end]/SecYear, τII_invariants[2:end]/1e6, color=:green, label="invariants, manually", linewidth=3)
lines!(ax, t_v[2:end]/SecYear, τII_invariants_2[2:end]/1e6, color=:yellow, label="invariants, elastic2", linewidth=3)

ax.xlabel= "Time (years)"
ax.ylabel= "Second invariant (MPa)"
display(fig)
