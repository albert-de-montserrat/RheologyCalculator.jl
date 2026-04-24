using RheologyCalculator, StaticArrays, ForwardDiff
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("../rheologies/RheologyDefinitions.jl")
include("tensor_helpers.jl")

function compute_stress_tensor(ε)
    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)
    plastic = DruckerPrager(1e6, 30, 0)

    # Maxwell visco-elasto-plastic model
    # elastic --- viscous --- plastic

    c  = SeriesModel(viscous, elastic, plastic)
    #c  = SeriesModel(viscous, elastic)
    
    # input variables (constant)
    εᵢⱼ    = Tuple(ε)
    τ0ᵢⱼ   = (zero_stress_tensor_2D(),)
    vars   = (; ε = εᵢⱼ, θ = 1.0e-20)
    # guess variables (we solve for these, differentiable)
    args   = (; τ = 0e0, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e8, P = 1.0e6, τ0 = τ0ᵢⱼ, P0 = (0.0,))

    x       = initial_guess_x(c, vars, args, others)
    char_τ  = plastic.C
    char_ε  = 1.0e-14
    xnorm   = normalisation_x(c, char_τ, char_ε)

    τII = solve(c, x, vars, others, verbose=true, xnorm0=xnorm)[1]
    τᵢⱼ = elastic_stress_history_2D(c, τII, vars.ε, τ0ᵢⱼ, others)[1]
    return SVector{length(τᵢⱼ)}(τᵢⱼ)
end

let 
    ε = @SVector [1e-14, -1e-14, 1e-16]
    D = ForwardDiff.jacobian(ε -> compute_stress_tensor(ε), ε)
end
