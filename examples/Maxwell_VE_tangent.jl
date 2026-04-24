using RheologyCalculator
using StaticArrays
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic

include("RheologyDefinitions.jl")

function compute_stress(ε)
    vars = (; ε = ε, )         # input variables (constant)
    args = (; τ = 2.0e3, )             # guess variables (we solve for these, differentiable)
    others = (; dt = 1.0e10, τ0 = (0e0, ), )       # other non-differentiable variables needed to evaluate the state functions

    x      = initial_guess_x(c, vars, args, others)
    char_τ = 1e6
    char_ε = vars.ε
    xnorm  = normalisation_x(c, char_τ, char_ε)

    x = solve(c, x, vars, others, verbose = false, xnorm=xnorm)[1] 
end

@inline compute_stress_tensor(ε::SVector{N, T}) where {N, T} = SVector{N, T}(compute_stress(ε[i]) for i in 1:N)

let 
    ε = @SVector rand(3)
    D = ForwardDiff.jacobian(ε -> compute_stress_tensor(ε), ε)
end
