using RheologyCalculator, ForwardDiff
using StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

function compute_stress_tensor(ε)
    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)
    c = SeriesModel(viscous, elastic)

    εᵢⱼ = Tuple(ε)
    τ0ᵢⱼ = (zero_stress_tensor_2D(),)
    vars = (; ε = εᵢⱼ, θ = 0.0)         # input variables (constant)
    args = (; τ = 2.0e3, P = 0.0)             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = τ0ᵢⱼ, P0 = (0.0,))       # other non-differentiable variables needed to evaluate the state functions

    x      = initial_guess_x(c, vars, args, others)
    char_τ = 1e6
    char_ε = 1.0e-14
    xnorm  = normalisation_x(c, char_τ, char_ε)

    τII = solve(c, x, vars, others, verbose=false, xnorm0=xnorm)[1]
    τᵢⱼ = elastic_stress_history_2D(c, τII, vars.ε, τ0ᵢⱼ, others)[1]
    return SVector{length(τᵢⱼ)}(τᵢⱼ)
end

let 
    ε = @SVector rand(3)
    D = ForwardDiff.jacobian(ε -> compute_stress_tensor(ε), ε)
end
