# This implements the Golchin yield function
using Test, LinearAlgebra
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using GLMakie 
using StaticArrays

Makie.inline!(true)
sc = (σ=1e7, t=1e10, L=1e3)

include("../rheologies/RheologyDefinitions.jl")
include("../rheologies/PorousDruckerPrager.jl")
include("tensor_helpers.jl")

# function stress_time(c, vars, x, xnorm, others; ntime = 200, dt = 1.0e8)
#     # Extract elastic stresses/pressure from solution vector
#     τ1      = zeros(ntime)
#     λ       = zeros(ntime)
#     P1      = zeros(ntime)
#     F1      = zeros(ntime)
    
#     mode2   = zeros(ntime)
#     t_v     = zeros(ntime)
#     τ_e     = (zero_stress_tensor_2D(),)
#     P_e     = (0.0e6,)
#     P1[1]   = P_e[1]
#     τ1[1]   = second_invariant_2D(τ_e[1])
#     x       = SA[τ1[1],0, P1[1]]
#     t       = 0.0
#     for i in 2:ntime
#         others = (; dt = dt, τ0 = τ_e, P0 = P_e)       # other non-differentiable variables needed to evaluate the state functions
        
#         x = RheologyCalculator.solve(c, x, vars, others, verbose = true, xnorm0=xnorm)
        
#         t += others.dt
        
#         τ_e = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
#         P_e = compute_pressure_elastic(c, x, others)
#         τ1[i] = second_invariant_2D(τ_e[1])
#         P1[i] = P_e[1]
#         F1[i] = compute_F(c.leafs[end], τ1[i], P1[i])   

#         t_v[i] = t
#     end

#     return t_v, τ1, P1, F1, mode2
# end

# The porosity rate is non-linear because bulk parameters (ηΦ) is a function of Φ
function porosity_rate(p̄, pf, λ̇, p̄0, pf0, KΦ, ηΦ, sinψ, Δt)
    dp̄dt    = (p̄ - p̄0) / Δt
    dpfdt   = (pf - pf0) / Δt
    dΦdt    = (dpfdt - dp̄dt)/KΦ + (pf - p̄)/ηΦ + λ̇*sinψ
    return dΦdt
end

# Equation of states solid and fluid
function EOS(p̄, pf, Φ, p̄0, pf0, Φ0, p)
    Ks, Kf, Δt = p.Ks, p.Kf, p.Δt
    dp̄dt   = (p̄ - p̄0) / Δt
    dpfdt   = (pf - pf0) / Δt
    dlnρfdt = dpfdt / Kf
    # Approximation in Yarushina ≈
    dPsdt   = ((p̄ - Φ*pf)/(1-Φ) - (p̄0 - Φ0*pf0)/(1-Φ0))/Δt
    # Exact, but non linear
    # dPsdt = dΦdt*(p̄ - pf*Φ)/(1-Φ)^2 + (dp̄dt - Φ*dpfdt - pf*dΦdt) / (1 - Φ)
    dlnρsdt = 1/Ks * dPsdt 
    return dlnρsdt, dlnρfdt
end

yield(τII, Pt, Pf, λ̇, c, sinϕ, cosϕ, ηvp) = τII - (Pt - Pf)*sinϕ - c*cosϕ - λ̇*ηvp

function residual_two_phase_trial(x, ε̇II_eff, divVs, divqD, p̄0, pf0, Φ0, c, Δt)

    τII, p̄, pf, λ̇ = x[1], x[2], x[3], x[4]
    Gs, Ks, KΦ, Kf = c.leafs[end-1].Gs, c.leafs[end-1].Ks , c.leafs[end-1].KΦ ,c.leafs[end-1].Kf  
    ηs, ηΦ = c.leafs[end-2].ηs, c.leafs[end-2].ηΦ
    C, sinϕ, cosϕ, sinψ, ηvp = c.leafs[end].C, c.leafs[end].sinϕ, c.leafs[end].cosϕ, c.leafs[end].sinΨ, c.leafs[end].η_vp
    eps   = -1e-13
    ηe    = Gs*Δt  
    ηve   = inv(1/ηs + 1/ηe)
    
    # Check yield
    f       = yield(τII, p̄, pf, λ̇, C, sinϕ, cosϕ, ηvp) 

    # Porosity rate
    dΦdt    = porosity_rate(p̄, pf, λ̇, p̄0, pf0, KΦ, ηΦ, sinψ, Δt)

    # Porosity implicitely updated
    Φ       = (Φ0 + dΦdt*Δt)

    # Form 3
    dpfdt   = (pf - pf0) / Δt
    dp̄dt    = (p̄ - p̄0) / Δt 
    dlnρfdt = dpfdt / Kf
    # # Exact, but non linear
    # dPsdt   = dΦdt*(p̄ - pf*Φ)/(1-Φ)^2 + (dp̄dt - Φ*dpfdt - pf*dΦdt) / (1 - Φ)
    # dlnρsdt = 1/Ks * dPsdt 
    # Approximation in Yarushina ≈
    dlnρsdt = 1/(1-Φ) *(dp̄dt - Φ*dpfdt) / Ks
    rp̄      = dlnρsdt - dΦdt/(1-Φ) + divVs
    rpf     = Φ*dlnρfdt + dΦdt     + Φ*divVs + divqD

    return @SVector[ 
        ε̇II_eff   -  τII/2/ηve - λ̇/2*(f>=eps),
        rp̄,
        rpf,
        (f - ηvp*λ̇)*(f>=eps) +  λ̇*(f<eps),
    ]
end

function stress_time(c, vars, x, xnorm, others; ntime = 30, dt = 1.0e10)

    # Extract elastic stresses/pressure from solution vector
    τ1      = zeros(ntime)
    λ̇1      = zeros(ntime)
    p̄1      = others.p̄0[1]*ones(ntime)
    p̄ᶠ1     = others.pᶠ0[1]*ones(ntime)
    Φ1      = others.Φ0[1]*ones(ntime)
    F1      = zeros(ntime)
    
    t_v     = zeros(ntime)
    τ_e     = (zero_stress_tensor_2D(),)
    p̄_e     = others.p̄0[1]
    pᶠ_e    = others.pᶠ0[1]
    Φ_e     = others.Φ0[1]
    p̄1[1]   = p̄_e[1]
    p̄ᶠ1[1]  = pᶠ_e[1]
    τ1[1]   = second_invariant_2D(τ_e[1])
    x       = SA[τ1[1],0, p̄1[1], p̄ᶠ1[1]]#, Φ1[1]]
    t       = 0.0
    for i in 2:ntime
        others = (; dt = dt, τ0 = τ_e, p̄0 = p̄_e, pᶠ0 = pᶠ_e, Φ0 = Φ_e)       # other non-differentiable variables needed to evaluate the state functions
        
        # x = RheologyCalculator.solve(c, x, vars, others, verbose = true, xnorm0=xnorm)

        ε̇     = [vars.ε[1], vars.ε[2], vars.ε[3]]  
        τ0    = [others.τ0[1][1], others.τ0[1][2], others.τ0[1][3]]
        Gs    = c.leafs[end-1].Gs
        ε̇_eff = ε̇ + τ0 / (2*Gs*dt)
        divVs = vars.θ
        divqD = vars.divqD
        p̄0, pf0, Φ0 = others.p̄0[1], others.pᶠ0[1], others.Φ0[1]

        ε̇II_eff = sqrt(0.5*(ε̇_eff[1]^2 + ε̇_eff[2]^2) + ε̇_eff[3]^2)

        @show "Step $(i)"
        nr0 = 1.0
        for iter=1:10
            r  = residual_two_phase_trial(x, ε̇II_eff, divVs, divqD, p̄0, pf0, Φ0, c, dt)
            nr = norm(r)
            if iter==1 nr0=nr end 
            @show nr, nr/nr0
            min(nr, nr/nr0)<1e-12 && break
            J = ForwardDiff.jacobian( x -> residual_two_phase_trial(x, ε̇II_eff, divVs, divqD, p̄0, pf0, Φ0, c, dt), x)
            x = x - J\r
        end  

        t += others.dt
        
        τ1[i] = x[1]
        p̄1[i] = x[2]
        p̄ᶠ1[i] = x[3]
        λ̇1[i] = x[4]
        # Φ1[i] = x[5]

        # In case porosity is post-processed (if eliminated from local iterations)
        Gs, Ks, KΦ, Kf = c.leafs[end-1].Gs, c.leafs[end-1].Ks , c.leafs[end-1].KΦ ,c.leafs[end-1].Kf  
        ηs, ηΦ = c.leafs[end-2].ηs, c.leafs[end-2].ηΦ
        C, sinϕ, cosϕ, sinψ, ηvp = c.leafs[end].C, c.leafs[end].sinϕ, c.leafs[end].cosϕ, c.leafs[end].sinΨ, c.leafs[end].η_vp
        dΦdt = porosity_rate(p̄1[i], p̄ᶠ1[i], λ̇1[i], p̄0, pf0, KΦ, ηΦ, sinψ, dt)[1]
        Φ1[i] = Φ0 + dΦdt*dt 
        
        τ_e  = (ε̇_eff.*τ1[i]./ε̇II_eff, )
        p̄_e  = (p̄1[i],)
        pᶠ_e = (p̄ᶠ1[i],)
        Φ_e  = (Φ1[i],)
        # τ_e   = elastic_stress_history_2D(c, x[1], vars.ε, τ_e, others)
        # P_e   = compute_pressure_elastic(c, x, others)
        # τ1[i] = second_invariant_2D(τ_e[1])
        # P1[i] = P_e[1]
        F1[i] = compute_F(c.leafs[end], τ1[i], p̄[i], pᶠ[i])   

        t_v[i] = t
    end

    return t_v, τ1, p̄1, p̄ᶠ1, Φ1, F1, λ̇1
end

c, x, xnorm, vars, args, others = let


    viscous = Poroviscosity(1e22/sc.σ/sc.t, 2e22/sc.σ/sc.t)
    elastic = Poroelasticity(3e10/sc.σ, 1e11/sc.σ, 1e10/sc.σ, 1e9/sc.σ)
    plastic = PorousDruckerPrager(; C=1e7/sc.σ, ϕ=35.0, ψ=10.0, η_vp=0e18/sc.σ/sc.t) 

    # Maxwell viscoelastic model
    c  = SeriesModel(viscous, elastic, plastic)

    # input variables (constant)
    @show "need a 2-phase dispatch of `vars_2D()` "
    # vars = vars_2D(0*7.0e-14, 7.0e-15)
    vars = (ε = (2.0e-15*sc.t, -2.0e-15*sc.t, 0.0), θ = 0.0, divqD = 0.0)
    @show vars

    # guess variables (we solve for these, differentiable)
    args = (; τ = 0.0e3, p̄ = 1e6/sc.σ, pᶠ = 1e6/sc.σ, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e5/sc.t, τ0 = (zero_stress_tensor_2D(),), p̄0 = (1e6/sc.σ), pᶠ0 = (1e6/sc.σ), Φ0 = 0.05)

    @show "need a 2-phase dispatch of `initial_guess_x()` "
    # x       = initial_guess_x(c, vars, args, others)
    x       = @SVector zeros(4)
    
    char_τ  = 1e7
    char_ε  = second_invariant_2D(vars.ε) + abs(vars.θ)
    xnorm   = normalisation_x(c, char_τ, char_ε)

    c, x, xnorm, vars, args, others
end

# Comute yield and potentials
τII = 0:1e6:1e8
p̄   = LinRange(-5e6, 2e8, 100)
pᶠ  = LinRange(-5e6, 1e8, 100)
F   = zeros(length(p̄), length(τII))
Q   = zeros(length(p̄), length(τII))
for i in CartesianIndices(F)
    F[i] = compute_F(c.leafs[end], τII[i[2]], p̄[i[1]], pᶠ[i[1]])
    Q[i] = compute_Q(c.leafs[end], τII[i[2]], p̄[i[1]], pᶠ[i[1]])
end

# Visualise
function figure()
    fig = Figure(fontsize = 20, size = (800, 800) )
    ax1 = Axis(fig[1,1], title="Dev. Stress",  xlabel=L"$t$ [yr]",  ylabel=L"$\tau_{II}$ [MPa]", xlabelsize=20, ylabelsize=20)
    ax2 = Axis(fig[2,1], title="Pressures",    xlabel=L"$t$ [yr]",  ylabel=L"$p^\text{eff}$ [MPa]", xlabelsize=20, ylabelsize=20)
    ax3 = Axis(fig[3,1], title="Plast. mult.", xlabel=L"$t$ [yr]",  ylabel=L"$\dot{\lambda}$ [1/s]", xlabelsize=20, ylabelsize=20)
    ax4 = Axis(fig[4,1], title="Porosity",     xlabel=L"$t$ [yr]",  ylabel=L"$\phi$ [-]",      xlabelsize=20, ylabelsize=20)

    SecYear = 3600 * 24 * 365.25
    t_v1, τ1, p̄1, pᶠ1, Φ1, F1, λ̇1  = stress_time(c, vars, x, xnorm, others; ntime = 30, dt = 1e10/sc.t)

    lines!(ax1, t_v1*sc.t / SecYear , τ1*sc.σ / 1.0e6,  color=:blue, label =  L"$\tau_{II}$")


    lines!(ax2, t_v1*sc.t / SecYear , p̄1*sc.σ  / 1.0e6, color=:red, label = L"$p^\text{eff}$")
    lines!(ax2, t_v1*sc.t / SecYear , pᶠ1*sc.σ / 1.0e6, color=:red, label = L"$p^\text{eff}$")

    lines!(ax3, t_v1*sc.t / SecYear , λ̇1/sc.t, color=:red, label = L"$p^\text{eff}$")

    lines!(ax4, t_v1*sc.t / SecYear , Φ1, color=:red, label = L"$p^\text{eff}$")

    # GLMakie.contour!(ax4, (p̄ .- pᶠ)/1e6, τII/1e6, F, levels = [0.001], color = :black)
    # GLMakie.contour!(ax4, (p̄ .- pᶠ)/1e6, τII/1e6, Q, levels = [0.001], color = :black, linestyle=:dash)
    # GLMakie.scatter!(ax4, P2/1e6, τ2/1e6, color = :red, label=L"1")
    # GLMakie.scatter!(ax4, P3/1e6, τ3/1e6, color = :blue, label=L"2")
    # GLMakie.scatter!(ax4, P1/1e6, τ1/1e6, color = :green, label=L"3")
    # axislegend(ax4, position=:lt)

    # GLMakie.save("./docs/assets/PorousDruckerPrager.png", fig)

    display(fig)
end

with_theme(figure, theme_latexfonts())
