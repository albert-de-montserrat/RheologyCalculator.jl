using ForwardDiff, StaticArrays
using GLMakie, MathTeXEngine
Makie.update_theme!( fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

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

viscous_strain_rate(η, τ)         = τ / (2 * η)
elastic_strain_rate(G, τ, τ0, dt) = (τ - τ0) / (2 * G * dt)
plastic_strain_rate(λ, τ, P, ψ)   = λ * ForwardDiff.derivative(τ -> compute_Q(τ, P, ψ), τ)

compute_Q(τ, P, ψ) = τ - P * sind(ψ)

function compute_F(τ, P, C, ϕ, λ)
    η_mult = 1.0  # Lagarange multiplier, value doesn't matter
    f      = τ - P * sind(ϕ) - C * cosd(ϕ) 
    F      = f*(f>=0) + η_mult*λ*(f<0)
    return F
end

function strain_rate(τ, τ0, dt, η, G, λ, P, ψ)
    ε_v = viscous_strain_rate(η, τ)
    ε_e = elastic_strain_rate(G, τ, τ0, dt)
    ε_p = plastic_strain_rate(λ, τ, P, ψ)
    ε = ε_v + ε_e + ε_p
    return ε
end

strain_rate_residual(ε, τ, τ0, dt, η, G, λ, P, ψ) = strain_rate(τ, τ0, dt, η, G, λ, P, ψ) - ε
F_residual(τ, P, C, ϕ, λ, G, dt, η) = compute_F(τ, P, C, ϕ, λ)

function residual_vector(x::SVector, ε, τ0, dt, η, G, P, ψ, C, ϕ)
    τ   = x[1]
    λ   = x[2]
    r_τ = strain_rate_residual(ε, τ, τ0, dt, η, G, λ, P, ψ)
    r_F = F_residual(τ, P, C, ϕ, λ, G, dt, η)
    return SA[r_τ, r_F]
end

function solve(ε, τ, τ0, dt, η, G, λ, P, ψ, C, ϕ; tol::Float64 = 1.0e-9, itermax = 1.0e4, verbose::Bool = false)

    it = 0
    er = Inf
    x  = SA[τ, λ]  # Initial guess
    α  = 1e0
    while er > tol #&& it<=1
        it += 1

        r = residual_vector(x, ε, τ0, dt, η, G, P, ψ, C, ϕ)  
        J = ForwardDiff.jacobian(x -> residual_vector(x, ε, τ0, dt, η, G, P, ψ, C, ϕ), x)
        Δx = J \ r
        α = 1 # bt_line_search(Δx, J, x, r, c, vars, others)
        x -= α .* Δx
        # check convergence
        er = mynorm(Δx, x)

        # @show r, x

        it > itermax && break
    end
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
    dt = 1e8
    ε  = 1e-14
    τ  = 0e6
    τ0 = 0
    η  = 1e22
    G  = 10e9
    λ  = 0
    P  = 1e6
    C  = 10e6
    ϕ  = 30
    ψ  = 0

    # Extract elastic stresses/pressure from solutio vector
    τv    = zeros(ntime)
    λv    = zeros(ntime)
    τv_an = zeros(ntime)
    tv    = zeros(ntime)
    t     = 0.0
    for i in 2:ntime
        sol          = solve(ε, τ, τ0, dt, η, G, λ, P, ψ, C, ϕ; verbose = true)
        τv[i], λv[i] = sol
        τ0           = sol[1]
        τ            = sol[1] # this is just a guess for the next iteration
        λ            = sol[2] # this is just a guess for the next iteration
        t           += dt
        τv_an[i]     = analytical_solution(ε, t, G, η, P, C, ϕ)
        tv[i]       = t
    end

    return tv, τv, τv_an
end


let
    tv, τv, τv_an = stress_time()

    fig = Figure(fontsize = 30, size = (800, 600) .* 2)
    ax1 = Axis(fig[1, 1], xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]", title=L"$$Stress - time")
    ax2 = Axis(fig[2, 1], xlabel = L"$t$ [kyr]", ylabel = L"$\tau$ [MPa]", title=L"$$Error")

    lines!(ax1, tv / SecYear / 1.0e3, τv_an / 1.0e6, color=:black, label = "analytical")
    scatter!(ax1, tv / SecYear / 1.0e3, τv / 1.0e6,  color=:red, label = "numerical")

    lines!(ax2, tv / SecYear / 1.0e3, (abs.(τv_an.-τv) ), color=:black)

    axislegend(ax1, position = :rb)
    display(fig)
end