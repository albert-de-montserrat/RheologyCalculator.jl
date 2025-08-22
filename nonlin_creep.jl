using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("examples/RheologyDefinitions.jl")

using GLMakie

function compute_viscosity(r, ε, T, P)
    (; n, A, E, V, R) = r
    η = A^(-1/n) * ε^((1-n)/n) * exp((E + P * V) / (n * R * T))
    return clamp(η, 1e16, 1e25)
end


function stress_strain(c, vars, x, others; dt = 1.0e8)
    ε = logrange(1e-20, 1e-14, 2_000)
    n = length(ε)
    # Extract elastic stresses/pressure from solutio vector
    τ1  = zeros(n)
    τ2  = zeros(n)
    η_eff1  = zeros(n)
    η_eff2  = zeros(n)
    for i in eachindex(ε)
        ## numeric
        vars   = (; ε = ε[i], )
        others = (; dt = dt, T = 1e3 + 273, P = 1.0e6, τ0 = 0e0)
        x      = solve(c, x, vars, others, verbose = true)
        η_eff1[i] = 0.5 * x[1] / ε[i]
        τ1[i]     = x[1]

        ## harmonic average
        η1    = compute_viscosity(c.leafs[1], ε[i], others.T, others.P)
        η2    = compute_viscosity(c.leafs[2], ε[i], others.T, others.P)
        # η_eff2[i] = 1/(1/η1 + 1/η2 + 1 / (c.leafs[3].G * dt))
        η_eff2[i] = 1/(1/η1 + 1/η2 )
        τ2[i] = 2 * η_eff2[i] * ε[i]
    end

    return ε, τ1, τ2, η_eff1, η_eff2
end

# Dry Olivine | Hirth & Kohlstedt (2003)
disl_data = (;
    n = 3.5,
    r = 0.0,
    A = 10^-15.96,
    E = 530.0e3,
    V = 13.0e-6,
    R = 8.314,
)
# Dry Olivine | Hirth & Kohlstedt (2003)
diff_data = (;
    n = 1e0,
    r = 0.0,
    p = -3,
    A = 10^-8.16,
    E = 375.0e3,
    V = 6e-6,
    R = 8.314,
)

dislocation = DislocationCreep(disl_data...)
diffusion   = DiffusionCreep(diff_data...)
elastic     = IncompressibleElasticity(74e9)
# Maxwell visco-elasto-plastic model
# c  = SeriesModel(dislocation, diffusion, elastic)
c  = SeriesModel(dislocation, diffusion)

# input variables (constant)
vars   = (; ε = 1.0e-14, )
# guess variables (we solve for these, differentiable)
args   = (; τ = 0e0, )
# other non-differentiable variables needed to evaluate the state functions
others = (; dt = 1.0e8, T = 1500e3 + 273, P = 3.88476e9, τ0 = 0e0)

# r=c.leafs[1]
# (; n, A, E, V, R) = r
# η = A^(-1/n) * ε[1]^((1-n)/n) * exp((E + 3.88476e9 * V) / (R * (1500 + 273)))
# η = A^(-1/n) * ε[end]^((1-n)/n) * exp((E + 3.88476e9 * V) / (R * (1500 + 273)))
# η = ε[end]^((1-n)/n) * A^(-1/n) *  exp((E + 3.88476e9 * V) / (R * (1500 + 273)))

x = initial_guess_x(c, vars, args, others)

ε, τ, τ2, η_eff1, η_eff2  = stress_strain(c, vars, x, others; dt = 50e3 * SecYear)

SecYear = 3600 * 24 * 365.25

fig = Figure(fontsize = 30, size = (800, 600) .* 2)
ax1  = Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = L"\dot{\varepsilon} [s^{-1}]", ylabel = L"\tau [MPa]")
ax2  = Axis(fig[2, 1], xscale = log10, yscale = log10, xlabel = L"\dot{\varepsilon} [s^{-1}]", ylabel = L"\eta_{\text{eff}}")

lines!(ax1, ε, τ  ./ 1e6,  color=:red, linewidth = 5, label = "numerical")
lines!(ax1, ε, τ2 ./ 1e6, color=:blue, linewidth = 5, label = "harmonic")

lines!(ax2, ε, η_eff1,  color=:red,  linewidth = 5, label = L"\eta_{\text{eff} \text{numerical}}")
lines!(ax2, ε, η_eff2,  color=:blue, linewidth = 5, label = L"\eta_{\text{eff} \text{harmonic}}")

axislegend(ax, position = :rb)
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [MPa]"
display(fig)


fig = Figure(fontsize = 30, size = (800, 600) .* 2)
ax1  = Axis(fig[1, 1], xscale = log10, xlabel = L"\dot{\varepsilon} [s^{-1}]", ylabel = L"\Delta\tau [MPa]")
ax2  = Axis(fig[2, 1], xscale = log10, xlabel = L"\dot{\varepsilon} [s^{-1}]", ylabel = L"\eta_{\text{numeric}_1} / \eta_{\text{harmonic}_2}")

lines!(ax1, ε, abs.(τ .- τ2) ./ 1e6, color=:black, linewidth = 5)
lines!(ax2, ε, η_eff1 ./ η_eff2, color=:black, linewidth = 5)

axislegend(ax, position = :rb)
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [MPa]"
display(fig)