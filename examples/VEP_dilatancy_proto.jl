using ForwardDiff, StaticArrays, LinearAlgebra
using GLMakie
#using MathTeXEngine
#Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

const SecYear = 3600 * 24 * 365.25

cancel_plastic_strain_rate = true

@generated function mynorm(x::SVector{N, T}, y::SVector{N, T}) where {N, T}
    return quote
        @inline
        v = zero(T)
        Base.@nexprs $N i -> begin
            xi = @inbounds x[i]
            yi = @inbounds y[i]
            v += !iszero(yi) * abs(xi / yi)
        end
        return v
    end
end

viscous_strain_rate(η, τ)                    = τ / (2 * η)
elastic_strain_rate(G, τ, τ0, dt)            = (τ - τ0) / (2 * G * dt)
volumetric_elastic_strain_rate(P, P0, K, dt) = (P - P0) / (K * dt)
function volumetric_plastic_strain_rate(λ, τ, P, C, ϕ, ψ, ηvp)   
    F   = compute_F(τ, P, C, ϕ, ηvp, λ)
    θ_p = λ   * ForwardDiff.derivative(P -> compute_Q(τ, P, ψ), P) 
    if cancel_plastic_strain_rate
        θ_p *= (F>-1e-8)
    end
    return θ_p
end

function plastic_strain_rate(λ, τ, P, C, ϕ, ψ,   ηvp)               
   F   = compute_F(τ, P, C, ϕ, ηvp, λ)
   ε_p = λ/2 * ForwardDiff.derivative(τ -> compute_Q(τ, P, ψ), τ)
   if cancel_plastic_strain_rate
       ε_p *= (F>-1e-8)
   end
   return ε_p
end

compute_Q(τ, P, ψ) = τ - P * sind(ψ)

function compute_F(τ, P, C, ϕ, ηvp, λ)
    η_mult = 1.0  # Lagarange multiplier, value doesn't matter
    f      = τ - P * sind(ϕ) - C * cosd(ϕ)
    f_vp   = τ - P * sind(ϕ) - C * cosd(ϕ) - ηvp*λ
    F      = f_vp
    return F
end

function strain_rate(τ, τ0, dt, η, G, λ, P, C, ϕ, ψ, ηvp)
    ε_v = viscous_strain_rate(η, τ)
    ε_e = elastic_strain_rate(G, τ, τ0, dt)
    ε_p = plastic_strain_rate(λ, τ, P, C, ϕ, ψ, ηvp)
    ε = ε_v + ε_e + ε_p
    return ε
end

function volumetric_strain_rate(τ, dt, λ, P, P0, K, C, ϕ, ψ, ηvp)
    θ_e = volumetric_elastic_strain_rate(P, P0, K, dt)
    θ_p = volumetric_plastic_strain_rate(λ, τ, P, C, ϕ, ψ, ηvp)
    θ   = θ_e + θ_p
    return θ
end
strain_rate_residual(ε, τ, τ0, dt, η, G, λ, P, C, ϕ, ψ, ηvp)         = strain_rate(τ, τ0, dt, η, G, λ, P, C, ϕ, ψ, ηvp) - ε
volumetric_strain_rate_residual(θ, τ, dt, λ, P, P0, K, C, ϕ, ψ, ηvp) = volumetric_strain_rate(τ, dt, λ, P, P0, K,  C, ϕ, ψ, ηvp) - θ
function F_residual(τ, P, C, ϕ, ηvp, λ, G, dt, η)                            
    F = compute_F(τ, P, C, ϕ, ηvp, λ) 
    return F*(F>-1e-8) - λ
end

function residual_vector(x::SVector, ε, τ0, dt, η, G, θ, P0, K, ψ, C, ϕ, ηvp)
    τ   = x[1]
    λ   = x[2]
    P   = x[3]
    r_F = F_residual(τ, P, C, ϕ, ηvp, λ, G, dt, η)
    r_τ = strain_rate_residual(ε, τ, τ0, dt, η, G, λ, P, C, ϕ, ψ, ηvp)
    r_θ = volumetric_strain_rate_residual(θ, τ, dt, λ, P, P0, K, C, ϕ, ψ, ηvp)
    return SA[r_τ, r_F, r_θ]
end

function solve(ε, τ, τ0, θ, dt, η, G, λ, P, P0, K, ψ, C, ϕ, ηvp; atol::Float64 = 1.0e-9, rtol::Float64 = 1.0e-9, itermax = 1.0e1, verbose::Bool = false)

    it = 0  
    er = Inf
    x  = SA[τ, λ, P]  # Initial guess
    α  = 1e0
    r   = residual_vector(x,ε, τ0, dt, η, G, θ, P0, K, ψ, C, ϕ, ηvp) 

    normalize_vec = SA[1.0, 1e3, 1.0]
    er0 = mynorm(r, normalize_vec)

    while er > atol && er/er0 > rtol #&& it<=1
        it += 1 

        J = ForwardDiff.jacobian(x -> residual_vector(x,ε, τ0, dt, η, G, θ, P0, K, ψ, C, ϕ, ηvp), x)
        
        #display(J)
        #error("stop")
        Δx = J \ r
        α = 1.0
        x -= α .* Δx
       
        # check convergence
        r   = residual_vector(x,ε, τ0, dt, η, G, θ, P0, K, ψ, C, ϕ, ηvp)  
        er  = mynorm(r, normalize_vec)

        it > itermax && break
      
        #=
        F = compute_F(x[1], x[3], C, ϕ, ηvp, x[2])

        println("               Iterations: $it, Error: $er, Error0: $er0, α = $α, F = $F, Error/Error0=$(er/er0)")
        println("                   r        = $r")
        println("                   x        = $x")
        println("                   Δx       = $Δx")
        println("                   Δx/(x+1) = $(Δx./ (x .+ 1.0))")
        println("                   r/(x+1)  = $(r./ (x .+ 1.0))")
        =#
    end
    #error("stop")
    if verbose
        println("Iterations: $it, Error: $er, α = $α")
    end
    return x
end

@inline function analytical_solution(ε, t, G, η, P, C, ϕ) 
    τ_ve = 2 * ε * η * (1 - exp(-G * t / η))
    τ_p  = P * sind(ϕ) + C * cosd(ϕ)
    return (τ_ve < τ_p) * τ_ve + (τ_ve >= τ_p) * τ_p
end

function stress_time()

    ntime = 2_000
    dt    = 1e8
    ε     = 1e-14
    θ     = 1e-18
    τ     = 1e3
    τ0    = 0
    η     = 1e22
    ηvp   = 1e20
    G     = 10e9
    K     = 30e9
    λ     = 0
    P     = 1e6
    P0    = P
    C     = 10e6
    ϕ     = 30
    ψ     = 30

    # Extract elastic stresses/pressure from solutio vector
    τv    = zeros(ntime)
    τv_pc = zeros(ntime)
    λv    = zeros(ntime)
    τv_an = zeros(ntime)
    Pv    = P0*ones(ntime)
    Pv_pc = P0*ones(ntime)
    tv    = zeros(ntime)
    t     = 0.0
    τ_pc  = 0.0
    p_pc  = P
    for i in 2:ntime
        sol                 = solve(ε, τ, τ0, θ, dt, η, G, λ, P, P0, K, ψ, C, ϕ, ηvp; verbose = true)
        τv[i], λv[i], Pv[i] = sol
        τ0                  = sol[1]
        τ                   = sol[1] # this is just a guess for the next iteration
        λ                   = sol[2] # this is just a guess for the next iteration
        P                   = sol[3] # this is just a guess for the next iteration
        P0                  = sol[3] 
        t                  += dt
        τv_an[i]            = analytical_solution(ε, t, G, η, P, C, ϕ)
        tv[i]               = t

        # Standard predictor-corrector
        η_ve       = (1/η + 1/(G*dt))^-1
        τ_pc0      = τ_pc
        p0         = p_pc
        τ_pc       = 2*η_ve * (ε + τ_pc0/(2*G*dt))
        f_p        = τ_pc - (C*cosd(ϕ) + p0*sind(ϕ))
        if f_p > 0
            λ    = f_p / (η_ve + ηvp + K*dt*sind(ϕ)*sind(ψ))
            εp   = λ*τ_pc/2/abs(τ_pc)
            divp = λ*sind(ψ)
            p_pc = p_pc + K*dt*divp
            f_c  = τ_pc - η_ve*λ - (C*cosd(ϕ) + p_pc*sind(ϕ) + λ*ηvp)   
            τ_pc = 2*η_ve * (ε + τ_pc0/(2*G*dt) - εp)
        end
        τv_pc[i]   = τ_pc
        Pv_pc[i]   = p_pc
    end

    return tv, τv, τv_pc, τv_an, Pv, Pv_pc
end


let
    tv, τv, τv_pc, τv_an, Pv, Pv_pc = stress_time()

    fig = Figure(fontsize = 30, size = (800, 600) .* 1)
    ax1 = Axis(fig[1, 1], xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]", title=L"$$Stress - time")
    ax2 = Axis(fig[2, 1], xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]", title=L"$$Error")
    ax3 = Axis(fig[3, 1], xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]", title=L"$$Pressure")

    step1 = 10
    step2 = 100
    lines!(ax1, tv / SecYear / 1.0e3, τv_an / 1.0e6, color=:black, label = "analytical")
    scatter!(ax1, tv[1:step1:end] / SecYear / 1.0e3, τv[1:step1:end]    / 1.0e6,  color=:red, label = "numerical")
    scatter!(ax1, tv[1:step2:end] / SecYear / 1.0e3, τv_pc[1:step2:end] / 1.0e6,  color=:blue, label = "predictor-corrector")

    lines!(ax2, tv / SecYear / 1.0e3, (abs.(τv_an.-τv) )/1e6, color=:black)
    lines!(ax2, tv / SecYear / 1.0e3, (abs.(τv_an.-τv_pc) )/1e6, color=:blue)

    scatter!(ax3, tv[1:step1:end] / SecYear / 1.0e3, Pv[1:step1:end]/1e6, color=:red)
    scatter!(ax3, tv[1:step2:end] / SecYear / 1.0e3, Pv_pc[1:step2:end]/1e6, color=:blue)

    axislegend(ax1, position = :rb)
    display(fig)
end