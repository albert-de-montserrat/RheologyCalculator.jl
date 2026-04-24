using RheologyCalculator, StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("RheologyDefinitions.jl")

function compute_stress(ε)
    viscous = LinearViscosity(1e22)
    elastic = IncompressibleElasticity(10e9)
    plastic = DruckerPrager(1e6, 30, 0)

    # Maxwell visco-elasto-plastic model
    # elastic --- viscous --- plastic

    c  = SeriesModel(viscous, elastic, plastic)
    #c  = SeriesModel(viscous, elastic)
    
    # input variables (constant)
    vars   = (; ε = ε, θ = 1.0e-20)
    # guess variables (we solve for these, differentiable)
    args   = (; τ = 0e0, λ = 0)
    # other non-differentiable variables needed to evaluate the state functions
    others = (; dt = 1.0e8, P = 1.0e6, τ0 = 0e0, P0 = 0.0)

    x       = initial_guess_x(c, vars, args, others)
    char_τ  = plastic.C
    char_ε  = vars.ε 
    xnorm   = normalisation_x(c, char_τ, char_ε)

    x = solve(c, x, vars, others, verbose = true, xnorm=xnorm)[1]
end

@inline compute_stress_tensor(ε::SVector{N, T}) where {N, T} = SVector{N, T}(compute_stress(ε[i]) for i in 1:N)

let 
    ε = @SVector [1e-14, -1e-14, 1e-16]
    D = ForwardDiff.jacobian(ε -> compute_stress_tensor(ε), ε)
end
