using TimerOutputs
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
        #α = bt_line_search_armijo(Δx, J, x, r, c, vars, others, α_min = 1.0e-8, c=0.9)
        α = bt_line_search(Δx, x, c, vars, others; α = 1.0, ρ = 0.5, lstol=0.9, α_min = 0.001) 
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


function solve_timed(c::AbstractCompositeModel, x::SVector, vars, others; tol::Float64 = 1.0e-9, itermax = 1.0e4, verbose::Bool = false)

    it = 0
    er = Inf
    local α
    to = TimerOutput()
    @timeit to "Newton-Raphson Iterations" while er > tol
        it += 1
        @timeit to "residual" r = compute_residual(c, x, vars, others)
        @timeit to "jacobian" J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
        @timeit to "backlash" Δx = J \ r
        @timeit to "line search" α = bt_line_search(Δx, J, x, r, c, vars, others)
        @timeit to "update solution" x -= α .* Δx
        # check convergence
        @timeit to "convergence check" er = mynorm(Δx, x)
        if verbose
            println("Iterations: $it, Error: $er, α = $α")
        end
        it > itermax && break
    end
    display(to)
    if verbose
        println("Iterations: $it, Error: $er, α = $α")
    end
    return x
end

function bt_line_search_armijo(Δx, J, x, r, composite, vars, others; α = 1.0, ρ = 0.5, c = 1.0e-4, α_min = 1.0e-8)

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


function bt_line_search(Δx, x, composite, vars, others; α = 1.0, ρ = 0.5, lstol=0.9, α_min = 1.0e-8)
    
    perturbed_x = @. x - α * Δx
    rnorm = compute_residual(composite, perturbed_x, vars, others)

    # Iterate unless step length becomes too small
    while α > α_min
        # Apply scaled update
        perturbed_x = @. x - α * Δx
        
        # Get updated residual
        rmnorm =   compute_residual(composite, perturbed_x, vars, others)

        # Check whether residual is sufficiently reduced
        if rmnorm <= lstol * rnorm
            break
        end

        # Bisect step length
        α *= ρ
    end

    return α
end

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
