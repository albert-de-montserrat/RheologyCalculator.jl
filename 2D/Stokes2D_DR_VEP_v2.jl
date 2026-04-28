# Initialisation
using Printf, Statistics, LinearAlgebra
using GLMakie
using TimerOutputs

using RheologyCalculator, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

import Statistics: mean

# second_invariant(ε) = √(0.5 * (ε[1]^2 + ε[2]^2) + ε[3]^2)
second_invariant(ε) = √(0.5 * (ε[1]^2 + ε[2]^2  + (-ε[1]-ε[2])^2) + ε[3]^2)

include("../rheologies/RheologyDefinitions.jl")

# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])
@views av4_harm(A) = 1.0./( 0.25.*(1.0./A[1:end-1,1:end-1].+1.0./A[2:end,1:end-1].+1.0./A[1:end-1,2:end].+1.0./A[2:end,2:end]) ) 

function dark_axis!(ax)
    ax.titlecolor = :white
    ax.xlabelcolor = :white
    ax.ylabelcolor = :white
    ax.xticklabelcolor = :white
    ax.yticklabelcolor = :white
    ax.xgridcolor = (:white, 0.18)
    ax.ygridcolor = (:white, 0.18)
    ax.xminorgridcolor = (:white, 0.08)
    ax.yminorgridcolor = (:white, 0.08)
    ax.leftspinecolor = :white
    ax.rightspinecolor = :white
    ax.topspinecolor = :white
    ax.bottomspinecolor = :white
    return ax
end

function dark_line_axis!(ax)
    dark_axis!(ax)
    ax.xgridcolor = (:white, 0.32)
    ax.ygridcolor = (:white, 0.32)
    ax.xgridwidth = 1.2
    ax.ygridwidth = 1.2
    return ax
end

function dark_colorbar!(cb)
    cb.labelcolor = :white
    cb.ticklabelcolor = :white
    cb.leftspinecolor = (:white, 0.25)
    cb.rightspinecolor = (:white, 0.25)
    cb.topspinecolor = (:white, 0.25)
    cb.bottomspinecolor = (:white, 0.25)
    return cb
end

# can be replaced by AD
function Gershgorin_Stokes2D_SchurComplement(ηc, ηv, γ, Δx, Δy, ncx  ,ncy)
        
    ηN    = ones(ncx-1, ncy)
    ηS    = ones(ncx-1, ncy)
    ηN[:,1:end-1] .= ηv[2:end-1,2:end-1]
    ηS[:,2:end-0] .= ηv[2:end-1,2:end-1]
    ηW    = ηc[1:end-1,:]
    ηE    = ηc[2:end-0,:]
    ebW   = γ[1:end-1,:] 
    ebE   = γ[2:end-0,:] 
    Cxx   = ones(ncx-1, ncy)
    Cxy   = ones(ncx-1, ncy)
    @. Cxx = abs.(ηN ./ Δy .^ 2) + abs.(ηS ./ Δy .^ 2) + abs.(ebE ./ Δx .^ 2 + (4 // 3) * ηE ./ Δx .^ 2) + abs.(ebW ./ Δx .^ 2 + (4 // 3) * ηW ./ Δx .^ 2) + abs.(-(-ηN ./ Δy - ηS ./ Δy) ./ Δy + (ebE ./ Δx + ebW ./ Δx) ./ Δx + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx)
    @. Cxy = abs.(ebE ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηN ./ (Δx .* Δy)) + abs.(ebE ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηS ./ (Δx .* Δy)) + abs.(ebW ./ (Δx .* Δy) + ηN ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy)) + abs.(ebW ./ (Δx .* Δy) + ηS ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy))
    
    Dx  = ones(ncx-1, ncy)
    @. Dx .= -(-ηN ./ Δy - ηS ./ Δy) ./ Δy + (ebE ./ Δx + ebW ./ Δx) ./ Δx + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx

    ηE    = ones(ncx, ncy-1)
    ηW    = ones(ncx, ncy-1)
    ηE[1:end-1,:] .= ηv[2:end-1,2:end-1]
    ηW[2:end-0,:] .= ηv[2:end-1,2:end-1]
    ηS    = ηc[:,1:end-1]
    ηN    = ηc[:,2:end-0]
    ebS  = γ[:,1:end-1] 
    ebN  = γ[:,2:end-0] 
    Cyy  = ones(ncx, ncy-1)
    Cyx  = ones(ncx, ncy-1)
    @. Cyy = abs.(ηE ./ Δx .^ 2) + abs.(ηW ./ Δx .^ 2) + abs.(ebN ./ Δy .^ 2 + (4 // 3) * ηN ./ Δy .^ 2) + abs.(ebS ./ Δy .^ 2 + (4 // 3) * ηS ./ Δy .^ 2) + abs.((ebN ./ Δy + ebS ./ Δy) ./ Δy + ((4 // 3) * ηN ./ Δy + (4 // 3) * ηS ./ Δy) ./ Δy - (-ηE ./ Δx - ηW ./ Δx) ./ Δx)
    @. Cyx = abs.(ebN ./ (Δx .* Δy) + ηE ./ (Δx .* Δy) - 2 // 3 * ηN ./ (Δx .* Δy)) + abs.(ebN ./ (Δx .* Δy) - 2 // 3 * ηN ./ (Δx .* Δy) + ηW ./ (Δx .* Δy)) + abs.(ebS ./ (Δx .* Δy) + ηE ./ (Δx .* Δy) - 2 // 3 * ηS ./ (Δx .* Δy)) + abs.(ebS ./ (Δx .* Δy) - 2 // 3 * ηS ./ (Δx .* Δy) + ηW ./ (Δx .* Δy))

    Dy  = ones(ncx, ncy-1)
    @. Dy .= (ebN ./ Δy + ebS ./ Δy) ./ Δy + ((4 // 3) * ηN ./ Δy + (4 // 3) * ηS ./ Δy) ./ Δy - (-ηE ./ Δx - ηW ./ Δx) ./ Δx

    λmaxVx = 1.0./Dx .* (Cxx .+ Cxy)
    λmaxVy = 1.0./Dy .* (Cyx .+ Cyy)

    return Dx, Dy, λmaxVx, λmaxVy
end

function solve_stress_RC(r, εxx, εyy, εxy, τ0xx, τ0yy, τ0xy, τII, P, dt)
    εᵢⱼ    = εxx, εyy, εxy
    τ0ᵢⱼ   = ((τ0xx, τ0yy, τ0xy), )
    vars   = (; ε = εᵢⱼ, θ = 0e0)                 # input variables (constant)
    args   = (; τ = τII, P = P)                   # guess variables (we solve for these, differentiable)
    others = (; dt = dt, τ0 = τ0ᵢⱼ, P0 = (0.0, )) # other non-differentiable variables needed to evaluate the state functions
    x      = initial_guess_x(r, vars, args, others)
    sol    = solve(r, x, vars, others, verbose = false)
    τII    = sol[1]
    
    # εII    = second_invariant(εᵢⱼ)
    # η      = τII / (2 * εII + eps())
    # τᵢⱼ    = @. 2 * εᵢⱼ * η

    ε_corr = effective_strain_rate_correction(r, εᵢⱼ, τ0ᵢⱼ, others)
    εeff   = εᵢⱼ .+ ε_corr
    εII    = second_invariant(εeff)
    η      = τII / (2 * εII + eps())
    # @show τII εII η
    τᵢⱼ    = @. 2 * εeff * η

    return τII, τᵢⱼ...
end

# @b solve_stress_RC($(r, Exxv[I], Eyyv[I], Exy[I], Txxv0[I], Tyyv0[I], Txy0[I], TIIv[I], Ptv[I], Δt)...)
# r = c[phases_v[I]]

# 2D Stokes routine
function Stokes2D_VEP(n)

    sc = (σ=1e6, t=1e10, L=1e3)

    # Maxwell viscoelastic model
    # elastic --- viscous

    viscous_matrix    = LinearViscosity(1e7)
    viscous_inclusion = LinearViscosity(1e-6)
    elastic           = IncompressibleElasticity(30e3)
    c_matrix          = SeriesModel(viscous_matrix, elastic)
    c_inclusion       = SeriesModel(viscous_inclusion, elastic)
    c                 = c_matrix, c_inclusion

    # Physics
    Lx, Ly   = 2e3/sc.L, 1e3/sc.L   # domain size
    radi     = 0.1e3/sc.L           # inclusion radius
    η0       = 1e23/sc.σ/sc.t       # viscous viscosity
    ηi       = 1e10/sc.σ/sc.t       # min/max inclusion viscosity
    G        = 3e10/sc.σ
    C        = 5e7/sc.σ
    Δt       = 4e10/sc.t
    εbg      =-1e-15*sc.t      # background strain-rate
    comp     = true                 
    K        = 5e10/sc.σ  
    ϕ        = 35.0 
    ψ        = 5.0   
    ηvp      = 2e20/sc.σ/sc.t    
    # Numerics
    ncx, ncy = 2*n*31, n*31   # numerical grid resolution
    nt       = 50           # time steps
    ϵ        = 1e-6         # tolerance
    iterMax  = 20000        # max number of iters
    nout     = 100           # check frequency
    c_fact   = 0.5          # damping factor
    dτ_local = true         # helps a little bit sometimes, sometimes not! 
    CFL_v    = 0.99         # CFL: can't make it larger
    γfact    = 20           # penalty: multiplier to the arithmetic mean of η
    rel_drop = 1e-1         # relative drop of velocity residual per PH iteration
    λ̇rel     = 1.075        # overrelaxation helps!
    # Preprocessing
    Δx, Δy  = Lx/ncx, Ly/ncy
    # Array initialisation
    Pt       = zeros(ncx  ,ncy  )
    Pt0      = zeros(ncx  ,ncy  ) 
    Ptv      = zeros(ncx+1,ncy+1)
    ΔPψ      = zeros(ncx  ,ncy  )
    ∇V       = zeros(ncx  ,ncy  )
    Vx       = zeros(ncx+1,ncy+2) 
    Vy       = zeros(ncx+2,ncy+1)
    dVx      = zeros(ncx-1,ncy  )
    dVy      = zeros(ncx  ,ncy-1)
    EIIc     = zeros(ncx  ,ncy  )
    EIIv     = zeros(ncx+1,ncy+1)
    Exx      = zeros(ncx  ,ncy  )
    Eyy      = zeros(ncx  ,ncy  )
    Exy      = zeros(ncx+1,ncy+1)
    Exxv     = zeros(ncx+1,ncy+1)
    Eyyv     = zeros(ncx+1,ncy+1)
    Exyc     = zeros(ncx  ,ncy  )
    TIIc     = zeros(ncx  ,ncy  )
    TIIv     = zeros(ncx+1,ncy+1)
    Txx      = zeros(ncx  ,ncy  )
    Tyy      = zeros(ncx  ,ncy  )
    Txy      = zeros(ncx+1,ncy+1)
    Txxv     = zeros(ncx+1,ncy+1)
    Tyyv     = zeros(ncx+1,ncy+1)
    Txyc     = zeros(ncx  ,ncy  )
    Txx0     = zeros(ncx  ,ncy  )
    Tyy0     = zeros(ncx  ,ncy  )
    Txy0     = zeros(ncx+1,ncy+1)
    Txxv0    = zeros(ncx+1,ncy+1)
    Tyyv0    = zeros(ncx+1,ncy+1)
    Txy0c    = zeros(ncx  ,ncy  )
    Fc       = zeros(ncx  ,ncy  ) 
    Fv       = zeros(ncx+1,ncy+1) 
    λ̇c       = zeros(ncx  ,ncy  )
    λ̇v       = zeros(ncx+1,ncy+1)
    λ̇_true_c = zeros(ncx  ,ncy  )
    λ̇_true_v = zeros(ncx+1,ncy+1)
    Rx       = zeros(ncx-1,ncy  )
    Ry       = zeros(ncx  ,ncy-1)
    Rp       = zeros(ncx  ,ncy  )
    Rx0      = zeros(ncx-1,ncy  )
    Ry0      = zeros(ncx  ,ncy-1)
    dVxdτ    = zeros(ncx-1,ncy  )
    dVydτ    = zeros(ncx  ,ncy-1)
    βVx      = zeros(ncx-1,ncy  )  # this disappears is dτ is not local
    βVy      = zeros(ncx  ,ncy-1)  # this disappears is dτ is not local
    cVx      = zeros(ncx-1,ncy  )  # this disappears is dτ is not local
    cVy      = zeros(ncx  ,ncy-1)  # this disappears is dτ is not local
    αVx      = zeros(ncx-1,ncy  )  # this disappears is dτ is not local
    αVy      = zeros(ncx  ,ncy-1)  # this disappears is dτ is not local
    ηb       = zeros(ncx  ,ncy  )
    ηc       = zeros(ncx  ,ncy  )
    ηv       = zeros(ncx+1,ncy+1)
    ηve_c    = zeros(ncx  ,ncy  )
    ηve_v    = zeros(ncx+1,ncy+1)
    ηvp_c    = zeros(ncx  ,ncy  )
    ηvp_v    = zeros(ncx+1,ncy+1)
    ηvep_c   = zeros(ncx  ,ncy  )
    ηvep_v   = zeros(ncx+1,ncy+1)
    ηc_sharp = zeros(ncx  ,ncy  )
    ηv_sharp = zeros(ncx+1,ncy+1)
    P_num    = zeros(ncx  ,ncy  )
    # Initialisation
    xce, yce = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx+2), LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, ncy+2)
    xc, yc   = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, ncx), LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, ncy)
    xv, yv   = LinRange(-Lx/2, Lx/2, ncx+1), LinRange(-Ly/2, Ly/2, ncy+1)
    # Multiple circles with various viscosities
    ηi    = (w=1/ηi, s=ηi) 
    x_inc = [0.0   0.2  -0.3 -0.4  0.0 -0.3 0.4  0.3  0.35 -0.1 ] 
    y_inc = [0.0   0.4   0.4 -0.3 -0.2  0.2 -0.2 -0.4 0.2  -0.4 ]
    r_inc = [radi  0.09  0.05 0.08 0.08  0.1 0.07 0.08 0.07 0.07] 
    η_inc = [ηi.s  ηi.w  ηi.w ηi.s ηi.w ηi.s ηi.w ηi.s ηi.s ηi.w]
    ηv_sharp   .= η0
    phases_c = [((x-x_inc[1])^2 + (y-y_inc[1]).^2 < r_inc[1]^2) ? 2 : 1   for x in xc, y in yc]
    phases_v = [((x-x_inc[1])^2 + (y-y_inc[1]).^2 < r_inc[1]^2) ? 2 : 1   for x in xv, y in yv]
    for inc in 1:1#eachindex(η_inc)
        ηv_sharp[(xv.-x_inc[inc]).^2 .+ (yv'.-y_inc[inc]).^2 .< r_inc[inc]^2 ] .= η_inc[inc]
    end
    ηc_sharp   .= η0
    for inc in 1:1#eachindex(η_inc)
        ηc_sharp[(xc.-x_inc[inc]).^2 .+ (yc'.-y_inc[inc]).^2 .< r_inc[inc]^2 ] .= η_inc[inc]
    end  
    # Harmonic averaging mimicking PIC interpolation
    ηc    .= av4_harm(ηv_sharp)
    ηv[2:end-1,2:end-1] .= av4_harm(ηc_sharp)
    ηv[1,:] .=  ηv[2,:]; ηv[end,:] .=  ηv[end-1,:]
    ηv[:,1] .=  ηv[:,2]; ηv[:,end] .=  ηv[:,end-1]
    # Effective viscosity
    ηve_c .= (1 ./ ηc .+ 1 ./ (G*Δt)).^-1
    ηve_v .= (1 ./ ηv .+ 1 ./ (G*Δt)).^-1
    # Bulk viscosity
    ηb    .= K .* Δt
    # Select γ
    γi   = γfact*mean(ηc).*ones(size(ηc))
    # (Pseudo-)compressibility
    γ_eff = zeros(size(ηb)) 
    if comp
        γ_num = γi.*ones(size(ηb)) * 1.0
        γ_phy = ηb
        γ_eff = ((γ_phy.*γ_num)./(γ_phy.+γ_num))
    else
        γ_eff .= γi
        γ_eff .= γ_eff
    end
    # Optimal pseudo-time steps - can be replaced by AD
    Dx, Dy, λmaxVx, λmaxVy = Gershgorin_Stokes2D_SchurComplement(ηve_c, ηve_v, γ_eff, Δx, Δy, ncx ,ncy)
    # Select dτ
    if dτ_local
        dτVx =  2.0./sqrt.(λmaxVx)*CFL_v
        dτVy =  2.0./sqrt.(λmaxVy)*CFL_v
    else
        dτVx =  2.0./sqrt.(maximum(λmaxVx))*CFL_v 
        dτVy =  2.0./sqrt.(maximum(λmaxVy))*CFL_v
    end
    βVx .= 2 .* dτVx ./ (2 .+ cVx.*dτVx)
    βVy .= 2 .* dτVy ./ (2 .+ cVy.*dτVy)
    αVx .= (2 .- cVx.*dτVx) ./ (2 .+ cVx.*dτVx)
    αVy .= (2 .- cVy.*dτVy) ./ (2 .+ cVy.*dτVy)
    # Initial condition
    Vx     .=   εbg.*xv .+    0*yce'
    Vy     .=     0*xce .- εbg.*yv'
    Vx[2:end-1,:] .= 0 # ensure non zero initial pressure residual
    Vy[:,2:end-1] .= 0 # ensure non zero initial pressure residual
    # Time
    Tii_evo = zeros(nt) 
    it_evo  = zeros(nt)
    TIIv2 = similar(TIIv)
    itg = 0
    err_evo_it, err_evo_V, err_evo_P = zeros(iterMax), zeros(iterMax), zeros(iterMax)
    to = TimerOutput()
    for it=1:nt
        Txx0 .= Txx; Tyy0 .= Tyy; Txy0 .= Txy; Txy0c .= Txyc;  Txxv0 .= Txxv; Tyyv0 .= Tyyv; Pt0 .= Pt
        # Iteration loop
        errVx0 = 1.0;  errVy0 = 1.0;  errPt0 = 1.0 
        errVx00= 1.0;  errVy00= 1.0; 
        iter=0; err=2*ϵ; err_evo_it .= 0.; err_evo_V .= 0.; err_evo_P .= 0.;
        @time for itPH = 1:100
            # Boundaries
            Vx[:,1] .= Vx[:,2]; Vx[:,end] .= Vx[:,end-1]
            Vy[1,:] .= Vy[2,:]; Vy[end,:] .= Vy[end-1,:]
            # Divergence
            ∇V    .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .+ (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy
            # Pressure on vertices
            Ptv[2:end-1,2:end-1] .= av(Pt)
            Ptv[1,:] .=  Ptv[2,:]; Ptv[end,:] .=  Ptv[end-1,:]
            Ptv[:,1] .=  Ptv[:,2]; Ptv[:,end] .=  Ptv[:,end-1]
            # Deviatoric strain rate
            Exx   .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .- 1.0/3.0.*∇V
            Eyy   .= (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy .- 1.0/3.0.*∇V
            Exy   .= 0.5.*((Vx[:,2:end] .- Vx[:,1:end-1])./Δy .+ (Vy[2:end,:] .- Vy[1:end-1,:])./Δx)
            Exxv[2:end-1,2:end-1] .= av(Exx)
            Eyyv[2:end-1,2:end-1] .= av(Eyy)
            Exyc  .= av(Exy)
            EIIc  .= sqrt.(0.5.*((Exx  .+ Txx0 ./(2*G*Δt)).^2 .+ (Eyy  .+ Tyy0 ./(2*G*Δt)).^2 .+ (.-(Exx  .+ Txx0 ./(2*G*Δt)).-(Eyy  .+ Tyy0 ./(2*G*Δt))).^2) .+ (Exyc .+ Txy0c./(2*G*Δt)).^2 )
            EIIv  .= sqrt.(0.5.*((Exxv .+ Txxv0./(2*G*Δt)).^2 .+ (Eyyv .+ Tyyv0./(2*G*Δt)).^2 .+ (.-(Exxv .+ Txxv0./(2*G*Δt)).-(Eyyv .+ Tyyv0./(2*G*Δt))).^2) .+ (Exy  .+ Txy0 ./(2*G*Δt)).^2 )
            # # Deviatoric stress
            # all(Txx   .== 2.0.*ηve_c.*(Exx  .+ Txx0 ./(2*G*Δt)))
            # all(Tyy   .== 2.0.*ηve_c.*(Eyy  .+ Tyy0 ./(2*G*Δt)))
            # all(Txy   .== 2.0.*ηve_v.*(Exy  .+ Txy0 ./(2*G*Δt)))
            # all(Txxv  .== 2.0.*ηve_v.*(Exxv .+ Txxv0./(2*G*Δt)))
            # all(Tyyv  .== 2.0.*ηve_v.*(Eyyv .+ Tyyv0./(2*G*Δt)))
            # all(Txyc  .== 2.0.*ηve_c.*(Exyc .+ Txy0c./(2*G*Δt)))
            # @timeit to "Arrays" begin
            #     Txx   .= 2.0.*ηve_c.*(Exx  .+ Txx0 ./(2*G*Δt))
            #     Tyy   .= 2.0.*ηve_c.*(Eyy  .+ Tyy0 ./(2*G*Δt))
            #     Txy   .= 2.0.*ηve_v.*(Exy  .+ Txy0 ./(2*G*Δt))
            #     Txxv  .= 2.0.*ηve_v.*(Exxv .+ Txxv0./(2*G*Δt))
            #     Tyyv  .= 2.0.*ηve_v.*(Eyyv .+ Tyyv0./(2*G*Δt))
            #     Txyc  .= 2.0.*ηve_c.*(Exyc .+ Txy0c./(2*G*Δt))
                
            #     TIIc  .= sqrt.(0.5.*(Txx.^2  .+ Tyy.^2  .+ (.-Txx.-Tyy).^2)   .+ Txyc.^2 )
            #     TIIv  .= sqrt.(0.5.*(Txxv.^2 .+ Tyyv.^2 .+ (.-Txxv.-Tyyv).^2) .+ Txy.^2 )
            # end
            # @timeit to "RC" begin
                # Threads.@threads 
            for I in CartesianIndices(Txxv)
                # centres
                if all(I.I .≤ size(TIIc))
                    TIIc[I], Txx[I], Tyy[I], Txyc[I] = solve_stress_RC(c[phases_c[I]], Exx[I], Eyy[I], Exyc[I], Txx0[I], Tyy0[I], Txy0c[I], TIIc[I], Pt[I], Δt)
                end
                # vertices
                TIIv[I], Txxv[I], Tyyv[I], Txy[I] = solve_stress_RC(c[phases_v[I]], Exxv[I], Eyyv[I], Exy[I], Txxv0[I], Tyyv0[I], Txy0[I], TIIv[I], Ptv[I], Δt)
            end
            # end
            # Plasticity
            λ̇c            .= 0.
            λ̇v            .= 0.
            Fc            .= TIIc .- C.*cosd(ϕ) .- Pt .*sind(ϕ)
            Fv            .= TIIv .- C.*cosd(ϕ) .- Ptv.*sind(ϕ)
            λ̇c[Fc.>0]     .= Fc[Fc.>0]./(ηve_c[Fc.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))      
            λ̇v[Fv.>0]     .= Fv[Fv.>0]./(ηve_v[Fv.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))      
            ηvep_c        .= ηve_c
            ηvep_v        .= ηve_v
            ηvp_c .= (TIIc.-λ̇c.*ηve_c) ./ (2 .* EIIc)
            ηvp_v .= (TIIv.-λ̇v.*ηve_v) ./ (2 .* EIIv)
            ηvep_c[Fc.>0] .= ηvp_c[Fc.>0]
            ηvep_v[Fv.>0] .= ηvp_v[Fv.>0]
            Txx   .= 2.0.*ηvep_c.*(Exx .+ Txx0./(2*G*Δt))
            Tyy   .= 2.0.*ηvep_c.*(Eyy .+ Tyy0./(2*G*Δt))
            Txy   .= 2.0.*ηvep_v.*(Exy .+ Txy0./(2*G*Δt))
            Txxv  .= 2.0.*ηvep_v.*(Exxv .+ Txxv0./(2*G*Δt))
            Tyyv  .= 2.0.*ηvep_v.*(Eyyv .+ Tyyv0./(2*G*Δt))
            Txyc  .= 2.0.*ηvep_c.*(Exyc .+ Txy0c./(2*G*Δt))
            ΔPψ   .= λ̇c.*sind(ψ).*K.*Δt
            # Check
            TIIc  .= sqrt.(0.5.*(Txx.^2  .+ Tyy.^2  .+ (.-Txx.-Tyy).^2)   .+ Txyc.^2 )
            TIIv  .= sqrt.(0.5.*(Txxv.^2 .+ Tyyv.^2 .+ (.-Txxv.-Tyyv).^2) .+ Txy.^2 )
            Fc    .= TIIc .- C.*cosd(ϕ) .- (Pt .+ λ̇c.*sind(ψ).*K.*Δt).*sind(ϕ)  .- ηvp.*λ̇c
            Fv    .= TIIv .- C.*cosd(ϕ) .- (Ptv.+ λ̇v.*sind(ψ).*K.*Δt).*sind(ϕ)  .- ηvp.*λ̇v
            # Residuals
            Rx    .= (.-(Pt[2:end,:] .- Pt[1:end-1,:])./Δx .- (ΔPψ[2:end,:] .- ΔPψ[1:end-1,:])./Δx .+ (Txx[2:end,:] .- Txx[1:end-1,:])./Δx .+ (Txy[2:end-1,2:end] .- Txy[2:end-1,1:end-1])./Δy)
            Ry    .= (.-(Pt[:,2:end] .- Pt[:,1:end-1])./Δy .- (ΔPψ[:,2:end] .- ΔPψ[:,1:end-1])./Δy .+ (Tyy[:,2:end] .- Tyy[:,1:end-1])./Δy .+ (Txy[2:end,2:end-1] .- Txy[1:end-1,2:end-1])./Δx)
            Rp    .= .-∇V .- comp*(Pt.-Pt0)./ηb 
            # Residual check
            errVx = norm(Rx); errVy = norm(Ry); errPt = norm(Rp)
            if itPH==1 errVx0=errVx; errVy0=errVy; errPt0=errPt; end
            err = maximum([min(errVx/errVx0, errVx), min(errVy/errVy0, errVy)]) #, min(errPt/errPt0, errPt)
            @printf("itPH = %02d iter = %06d iter/nx = %03d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, Rp=%1.3e] - max(F) = %1.2e \n", itPH, iter, iter/ncx, err, min(errVx/errVx0, errVx), min(errVy/errVy0, errVy), min(errPt/errPt0, errPt), max(maximum(Fc), maximum(Fv)))
            if (err<ϵ) break end
            # Set tolerance of velocity solve proportional to residual
            ϵ_vel = err*rel_drop
            itPT  = 0.
            while (err>ϵ_vel && itPT<=iterMax)
                iter   += 1 
                itPT   += 1
                itg    += 1
                # Pseudo-old dudes 
                Rx0   .= Rx
                Ry0   .= Ry
                # Boundaries
                Vx[:,1] .= Vx[:,2]; Vx[:,end] .= Vx[:,end-1]
                Vy[1,:] .= Vy[2,:]; Vy[end,:] .= Vy[end-1,:]
                # Divergence 
                ∇V    .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .+ (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy
                Rp    .= .-∇V .- comp*(Pt.-Pt0)./ηb 
                # Deviatoric strain rate
                Exx   .= (Vx[2:end,2:end-1] .- Vx[1:end-1,2:end-1])./Δx .- 1.0/3.0.*∇V
                Eyy   .= (Vy[2:end-1,2:end] .- Vy[2:end-1,1:end-1])./Δy .- 1.0/3.0.*∇V
                Exy   .= 0.5.*((Vx[:,2:end] .- Vx[:,1:end-1])./Δy .+ (Vy[2:end,:] .- Vy[1:end-1,:])./Δx)
                Exxv[2:end-1,2:end-1] .= av(Exx)
                Eyyv[2:end-1,2:end-1] .= av(Eyy)
                Exyc  .= av(Exy)
                EIIc  .= sqrt.(0.5.*((Exx  .+ Txx0 ./(2*G*Δt)).^2 .+ (Eyy  .+ Tyy0 ./(2*G*Δt)).^2 .+ (.-(Exx  .+ Txx0 ./(2*G*Δt)).-(Eyy  .+ Tyy0 ./(2*G*Δt))).^2) .+ (Exyc .+ Txy0c./(2*G*Δt)).^2 )
                EIIv  .= sqrt.(0.5.*((Exxv .+ Txxv0./(2*G*Δt)).^2 .+ (Eyyv .+ Tyyv0./(2*G*Δt)).^2 .+ (.-(Exxv .+ Txxv0./(2*G*Δt)).-(Eyyv .+ Tyyv0./(2*G*Δt))).^2) .+ (Exy  .+ Txy0 ./(2*G*Δt)).^2 )
                # # Deviatoric stress
                # Txx   .= 2.0.*ηve_c.*(Exx  .+ Txx0 ./(2*G*Δt)) 
                # Tyy   .= 2.0.*ηve_c.*(Eyy  .+ Tyy0 ./(2*G*Δt)) 
                # Txy   .= 2.0.*ηve_v.*(Exy  .+ Txy0 ./(2*G*Δt))
                # Txxv  .= 2.0.*ηve_v.*(Exxv .+ Txxv0./(2*G*Δt))
                # Tyyv  .= 2.0.*ηve_v.*(Eyyv .+ Tyyv0./(2*G*Δt))
                # Txyc  .= 2.0.*ηve_c.*(Exyc .+ Txy0c./(2*G*Δt))
                # TIIc  .= sqrt.(0.5.*(Txx.^2  .+ Tyy.^2  .+ (.-Txx.-Tyy).^2)   .+ Txyc.^2 )
                # TIIv  .= sqrt.(0.5.*(Txxv.^2 .+ Tyyv.^2 .+ (.-Txxv.-Tyyv).^2) .+ Txy.^2 )
                Threads.@threads for I in CartesianIndices(Txxv)
                    # centres
                    if all(I.I .≤ size(TIIc))
                        TIIc[I], Txx[I], Tyy[I], Txyc[I] = solve_stress_RC(c[phases_c[I]], Exx[I], Eyy[I], Exyc[I], Txx0[I], Tyy0[I], Txy0c[I], TIIc[I], Pt[I], Δt)
                    end
                    # vertices
                    TIIv[I], Txxv[I], Tyyv[I], Txy[I] = solve_stress_RC(c[phases_v[I]], Exxv[I], Eyyv[I], Exy[I], Txxv0[I], Tyyv0[I], Txy0[I], TIIv[I], Ptv[I], Δt)
                end
                # Plasticity
                Fc              .= TIIc .- C.*cosd(ϕ) .- Pt .*sind(ϕ)
                Fv              .= TIIv .- C.*cosd(ϕ) .- Ptv.*sind(ϕ)
                λ̇_true_c        .= 0.
                λ̇_true_v        .= 0.
                λ̇_true_c[Fc.>0] .= Fc[Fc.>0]./(ηve_c[Fc.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))      
                λ̇_true_v[Fv.>0] .= Fv[Fv.>0]./(ηve_v[Fv.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))
                λ̇c              .= λ̇rel*λ̇_true_c .+ (1-λ̇rel).*λ̇c
                λ̇v              .= λ̇rel*λ̇_true_v .+ (1-λ̇rel).*λ̇v 
                ηvep_c .= ηve_c
                ηvep_v .= ηve_v
                ηvp_c  .= (TIIc.-λ̇c.*ηve_c) ./ (2 .* EIIc)
                ηvp_v  .= (TIIv.-λ̇v.*ηve_v) ./ (2 .* EIIv)
                ηvep_c[Fc.>0] .= ηvp_c[Fc.>0]
                ηvep_v[Fv.>0] .= ηvp_v[Fv.>0] 
                Txx    .= 2.0.*ηvep_c.*(Exx  .+ Txx0 ./(2*G*Δt)) 
                Tyy    .= 2.0.*ηvep_c.*(Eyy  .+ Tyy0 ./(2*G*Δt)) 
                Txy    .= 2.0.*ηvep_v.*(Exy  .+ Txy0 ./(2*G*Δt))
                ΔPψ    .= λ̇c.*sind(ψ).*K.*Δt
                # Residuals
                P_num  .= γ_eff .* Rp
                Rx     .= (1.0./Dx).*(.-(P_num[2:end,:] .- P_num[1:end-1,:])./Δx .- (Pt[2:end,:] .- Pt[1:end-1,:])./Δx .- (ΔPψ[2:end,:] .- ΔPψ[1:end-1,:])./Δx .+ (Txx[2:end,:] .- Txx[1:end-1,:])./Δx .+ (Txy[2:end-1,2:end] .- Txy[2:end-1,1:end-1])./Δy)
                Ry     .= (1.0./Dy).*(.-(P_num[:,2:end] .- P_num[:,1:end-1])./Δy .- (Pt[:,2:end] .- Pt[:,1:end-1])./Δy .- (ΔPψ[:,2:end] .- ΔPψ[:,1:end-1])./Δy .+ (Tyy[:,2:end] .- Tyy[:,1:end-1])./Δy .+ (Txy[2:end,2:end-1] .- Txy[1:end-1,2:end-1])./Δx)
                # Damping-pong
                dVxdτ  .= αVx.*dVxdτ .+ Rx
                dVydτ  .= αVy.*dVydτ .+ Ry
                # PT updates
                Vx[2:end-1,2:end-1] .+= dVxdτ.*βVx.*dτVx 
                Vy[2:end-1,2:end-1] .+= dVydτ.*βVy.*dτVy 
                # Residual check
                if mod(iter, nout)==0
                    errVx = norm(Dx.*Rx); errVy = norm(Dy.*Ry)
                    if iter==nout errVx00=errVx; errVy00=errVy; end
                    err = maximum([errVx./errVx00, errVy./errVy00])
                    err_evo_V[iter] = errVx/errVx00; err_evo_P[iter] = errPt/errPt0; err_evo_it[iter] =  iter
                    dVx .= dVxdτ.*βVx.*dτVx
                    dVy .= dVydτ.*βVy.*dτVy
                    @printf("it = %d, iter = %d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e] \n", it, iter, err, errVx/errVx00, errVy/errVy00)
                    # λminV  = abs.((sum(dVx.*(Rx .- Rx0))) + abs.((sum(dVy.*(Ry .- Ry0))) )/ ( sum(dVx.*dVx)) + sum(dVy.*dVy) ) 
                    λminV  = abs(  sum(dVx.*(Rx .- Rx0)) + sum(dVy.*(Ry .- Ry0))  ) / (sum(dVx.*dVx) .+ sum(dVy.*dVy))
                    cVx .= 2*sqrt.(λminV)*c_fact
                    cVy .= 2*sqrt.(λminV)*c_fact
                    # Optimal pseudo-time steps - can be replaced by AD
                    Dx, Dy, λmaxVx, λmaxVy = Gershgorin_Stokes2D_SchurComplement(ηve_c, ηve_v, γ_eff, Δx, Δy, ncx ,ncy)
                    # Select dτ
                    if dτ_local
                        dτVx =  2.0./sqrt.(λmaxVx)*CFL_v
                        dτVy =  2.0./sqrt.(λmaxVy)*CFL_v
                    else
                        dτVx =  2.0./sqrt.(maximum(λmaxVx))*CFL_v 
                        dτVy =  2.0./sqrt.(maximum(λmaxVy))*CFL_v
                    end
                    βVx .= 2 .* dτVx ./ (2 .+ cVx.*dτVx)
                    βVy .= 2 .* dτVy ./ (2 .+ cVy.*dτVy)
                    αVx .= (2 .- cVx.*dτVx) ./ (2 .+ cVx.*dτVx)
                    αVy .= (2 .- cVy.*dτVy) ./ (2 .+ cVy.*dτVy)
                end
            end
            Pt .+= γ_eff.*Rp
        end

        Tii_evo[it] = maximum(TIIc)
        it_evo[it]  = iter/ncx

        # Plotting
        EIIc  .= sqrt.(0.5.*((Exx).^2 .+ (Eyy).^2 .+ (.-(Exx).-(Eyy)).^2) .+ (Exyc).^2 )
        fig = Figure(size = (1100, 850), backgroundcolor = :black)

        ax1 = Axis(fig[1, 1], title = "EII", xlabel = "x", ylabel = "y", aspect = DataAspect(), backgroundcolor = :black)
        dark_axis!(ax1)
        hm1 = heatmap!(ax1, xc, yc, log10.(EIIc ./ sc.t), colormap = :inferno)
        xlims!(ax1, -Lx / 2, Lx / 2)
        dark_colorbar!(Colorbar(fig[1, 2], hm1))

        ax2 = Axis(fig[1, 3], title = "Pt", xlabel = "x", ylabel = "y", aspect = DataAspect(), backgroundcolor = :black)
        dark_axis!(ax2)
        hm2 = heatmap!(ax2, xc, yc, Pt .* sc.σ, colormap = :inferno)
        xlims!(ax2, -Lx / 2, Lx / 2)
        dark_colorbar!(Colorbar(fig[1, 4], hm2))

        ax3 = Axis(fig[2, 1:2], xlabel = "time", ylabel = "mean dev. stress", backgroundcolor = :black)
        dark_line_axis!(ax3)
        lines!(ax3, 1:it, Tii_evo[1:it] .* sc.σ, color = :deepskyblue2, linewidth = 3)

        ax4 = Axis(fig[2, 3], title = "ηc", xlabel = "x", ylabel = "y", aspect = DataAspect(), backgroundcolor = :black)
        dark_axis!(ax4)
        hm4 = heatmap!(ax4, xc, yc, log10.(ηc), colormap = :inferno)
        xlims!(ax4, -Lx / 2, Lx / 2)
        dark_colorbar!(Colorbar(fig[2, 4], hm4))

        colgap!(fig.layout, 35)

        display(fig)
        @show iter/ncx
        @show itg

    end
    n   = length(ηc)
    @show η_h = 1.0 / sum(1.0/n ./ηc)
    @show η_g = exp( sum( 1.0/n*log.(ηc)))
    @show η_a = mean(ηc)
    # @show to
    return
end

Stokes2D_VEP(2)


#### V-E case

# ────────────────────────────────────────────────────────────────────
#                            Time                    Allocations      
#                   ───────────────────────   ────────────────────────
# Tot / % measured:      8.07s /   0.8%           30.2GiB /   0.0%    

# Section   ncalls     time    %tot     avg     alloc    %tot      avg
# ────────────────────────────────────────────────────────────────────
# RC            58   63.7ms   96.2%  1.10ms     0.00B     - %    0.00B
# Arrays        58   2.55ms    3.8%  44.0μs     0.00B     - %    0.00B
# ────────────────────────────────────────────────────────────────────