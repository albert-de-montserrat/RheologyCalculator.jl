"""
    x = solve(c::AbstractCompositeModel, x::SVector, vars, others; tol = 1.0e-9, itermax = 1e4, verbose=true)

Solve the system of equations defined by the composite model `c` using a Newton-Raphson method.
"""
function solve(c::AbstractCompositeModel, x::SVector, vars, others; tol::Float64 = 1.0e-9, itermax = 1.0e4, verbose::Bool = false)

    it = 0
    er = Inf
    local α
    while er > tol
        it += 1
        r = compute_residual(c, x, vars, others)
        J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
        Δx = J \ r
        α = bt_line_search(Δx, J, x, r, c, vars, others)
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

function bt_line_search(Δx, J, x, r, composite, vars, others; α = 1.0, ρ = 0.5, c = 1.0e-4, α_min = 1.0e-8)

    perturbed_x = @. x - α * Δx
    perturbed_r = compute_residual(composite, perturbed_x, vars, others)

    J_times_Δx = J * Δx
    while sqrt(sum(perturbed_r .^ 2)) > sqrt(sum((r + (c * α * (J_times_Δx))) .^ 2))
        α *= ρ
        if α < α_min
            α = α_min
            break
        end
        perturbed_x = @. x - α * Δx
        perturbed_r = compute_residual(composite, perturbed_x, vars, others)
    end
    return α
end

@generated function mynorm(x::SVector{N,T}, y::SVector{N,T}) where {N, T}
    quote 
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
