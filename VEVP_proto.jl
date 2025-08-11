using ForwardDiff, StaticArrays
using GLMakie

const SecYear = 3600 * 24 * 365.25

@generated function mynorm(x::SVector{N, T}, y::SVector{N, T}) where {N, T}
    return quote
        @inline
        v = zero(T)
        Base.@nexprs $N i -> begin
            xi = @inbounds x[i]
            yi = @inbounds y[i]
            v += !iszero(xi) * abs(xi / yi)
        end
        return v
    end
end

viscous_stress(η, ε)              = 2 * η * ε
viscous_strain_rate(η, τ)         = τ / (2 * η)
elastic_strain_rate(G, τ, τ0, dt) = (τ - τ0) / (2 * G * dt)
plastic_strain_rate(λ, τ, P, ψ)   = λ / 2 * ForwardDiff.derivative(τ -> compute_Q(τ, P, ψ), τ)

compute_Q(τ, P, ψ) = τ - P * sind(ψ)

function compute_F(τ, P, C, ϕ, λ, ηve)
    F = τ - P * sind(ϕ) - C * cosd(ϕ) #- λ * ηve
    F *= (F > 0)
    return F - λ * 1
end

function strain_rate(τ, τ0, dt, η, G, ε_p)
    ε_v = viscous_strain_rate(η, τ)
    ε_e = elastic_strain_rate(G, τ, τ0, dt)
    ε = ε_v + ε_e + ε_p
    return ε
end

strain_rate_residual(ε, τ, τ0, dt, η, G, ε_p) = strain_rate(τ, τ0, dt, η, G, ε_p) - ε
F_residual(τ, P, C, ϕ, λ, ηve)                = compute_F(τ, P, C, ϕ, λ, ηve)
parallel_strain_rate_residual(τ, τpl, ε, η)   = τpl + viscous_stress(η, ε) - τ
plastic_stress_residual(τ, P, ε, λ, ψ)        = plastic_strain_rate(λ, τ, P, ψ) - ε


function residual_vector(x::SVector, ε, τ0, dt, η, ηvp, G, P, ψ, C, ϕ)
    τ     = x[1]
    εp    = x[2]
    τpl   = x[3]
    λ     = x[4]
    ηve   = 1/(1/(η) + 1/(G*dt))
    r_τ   = strain_rate_residual(ε, τ, τ0, dt, η, G, εp)
    r_εpl = parallel_strain_rate_residual(τ, τpl, εp, ηvp)
    r_τpl = plastic_stress_residual(τpl, P, εp, λ, ψ)
    r_F   = F_residual(τpl, P, C, ϕ, λ, ηve)
    return SA[r_τ, r_εpl, r_τpl, r_F]
end

function solver(ε, τ, τ0, dt, η, ηvp, G, λ, P, ψ, C, ϕ; tol::Float64 = 1.0e-9, itermax = 1.0e4, verbose::Bool = false)

    it = 0
    er = Inf
    x  = SA[τ, ε, 0, 0]  # Initial guess
    α  = 1e0
    while er > tol
        it += 1

        r = residual_vector(x, ε, τ0, dt, η, ηvp, G, P, ψ, C, ϕ)
        J = ForwardDiff.jacobian(x -> residual_vector(x, ε, τ0, dt, η, ηvp, G, P, ψ, C, ϕ), x)
        Δx = J \ r
        α  = 1 # bt_line_search(Δx, J, x, r, c, vars, others)
        x -= α .* Δx
        # check convergence
        er = mynorm(Δx, x)

        it > itermax && break
    end
    if verbose
        println("Iterations: $it, Error: $er, α = $α")
    end
    return x
end

@inline analytical_solution(ε, t, G, η) = 2 * ε * η * (1 - exp(-G * t / η))

function stress_time()

    ntime = 1_500
    dt    = 1e8
    ε     = 1e-14
    τ     = 2e3
    τ0    = 0
    η     = 1e22
    ηvp   = 1e20
    G     = 10e9
    λ     = 1e-12
    P     = 1e6
    C     = 10e6
    ϕ     = 30
    ψ     = 0

    # Extract elastic stresses/pressure from solutio vector
    τv    = zeros(ntime)
    τv_pc = zeros(ntime)
    λv    = zeros(ntime)
    τv_an = zeros(ntime)
    tv    = zeros(ntime)
    t     = 0.0
    τ_pc  = 0.0
    sol_t = Any[]
    for i in 2:ntime
        sol       = solver(ε, τ, τ0, dt, η, ηvp, G, λ, P, ψ, C, ϕ; verbose = false)
        push!(sol_t, sol)
        τv[i]     = sol[1]
        τ         = sol[1]
        τ0        = sol[1]
        t        += dt
        τv_an[i]  = analytical_solution(ε, t, G, η)
        tv[i]     = t


        # Standard predictor-corrector
        η_ve       = (1/η + 1/(G*dt))^-1
        τ_pc0      = τ_pc
        τ_pc       = 2*η_ve * (ε + τ_pc0/(2*G*dt))
        f_p        = τ_pc - (C*cosd(ϕ) + P*sind(ϕ))
        if f_p > 0
            λ    = f_p / (η_ve + ηvp)
            εp   = λ*τ_pc/2/abs(τ_pc)
            f_c  = τ_pc - η_ve*λ - (C*cosd(ϕ) + P*sind(ϕ) + λ*ηvp)   
            @show f_c
            τ_pc = 2*η_ve * (ε + τ_pc0/(2*G*dt) - εp)
        end
        τv_pc[i]   = τ_pc

    end

    return tv, τv, τv_pc, τv_an, sol_t   
end


tv, τv, τv_pc, τv_an, sol_t_proto = stress_time()

fig = Figure(fontsize = 30, size = (800, 600) .* 2)
ax  = Axis(fig[1, 1], xlabel = "t [kyr]", ylabel = L"\tau [MPa]")
# ax2 = Axis(fig[2, 1], xlabel = "t [kyr]", ylabel = L"\tau [MPa]")

step1 = 10
step2 = 100
lines!(ax, tv / SecYear / 1.0e3, τv_an / 1.0e6, color=:black, label = "analytical")
scatter!(ax, tv[1:step1:end] / SecYear / 1.0e3, τv[1:step1:end]    / 1.0e6,  color=:red, label = "numerical")
scatter!(ax, tv[1:step2:end] / SecYear / 1.0e3, τv_pc[1:step2:end] / 1.0e6,  color=:blue, label = "predictor-corrector")

# lines!(ax2, tv / SecYear / 1.0e3, log10.(abs.(τ_an.-τ) ./ τ_an), color=:black)

# axislegend(ax, position = :rb)
#title!(ax,"Burgers model")
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [MPa]"

axislegend(ax, position = :rb)

display(fig)