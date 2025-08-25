using TimerOutputs
"""
    x = solve(c::AbstractCompositeModel, x::SVector, vars, others; atol = 1.0e-9, rtol = 1.0e-9, itermax = 1e4, verbose=true)

Solve the system of equations defined by the composite model `c` using a Newton-Raphson method.
Optional parameters:
- `xnorm`: Normalization vector for the solution vector `x`.
- `atol`: Absolute tolerance for convergence (default: 1.0e-9)
- `rtol`: Relative tolerance for convergence (default: 1.0e-9)
- `itermax`: Maximum number of iterations (default: 1e4)
- `verbose`: If true, print convergence information (default: false)
"""
function solve(c::AbstractCompositeModel, x::SVector, vars, others; xnorm=nothing, atol::Float64 = 1.0e-9, rtol::Float64 = 1.0e-9, itermax = 1.0e4, verbose::Bool = false)

    if isnothing(xnorm)
        xnorm = x.*0 .+ 1     # set normalization vector to 1.0
    end
    r  = compute_residual(c, x, vars, others)   # initial residual

    it  = 0
    er  = Inf
    er0 = mynorm(r, xnorm)

    local α
    while er > atol && er > rtol * er0 
        it += 1
        
        J = ForwardDiff.jacobian(y -> compute_residual(c, y, vars, others), x)
        Δx = J \ -r
        #α = bt_line_search_armijo(Δx, J, x, xnorm, c, vars, others, α_min = 1.0e-8, c=0.9)
        α = bt_line_search(Δx, x, c, vars, others, xnorm; α = 1.0, ρ = 0.5, lstol=0.95, α_min = 0.1) 
        x += α .* Δx

        # check convergence
        r  = compute_residual(c, x, vars, others)
        er = mynorm(r, xnorm)

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

function bt_line_search_armijo(Δx, J, x, xnorm, composite, vars, others; α = 1.0, ρ = 0.5, c = 1.0e-4, α_min = 1.0e-8)

    r           = compute_residual(composite, x, vars, others)
    rnorm       = mynorm(r, xnorm)
    J_times_Δx  = J * Δx
    while α > α_min
        perturbed_x = @. x + α * Δx

        perturbed_r     = compute_residual(composite, perturbed_x, vars, others)
        perturbed_rnorm = mynorm(perturbed_r, xnorm)  

        armijo_condition = perturbed_rnorm^2 ≤ rnorm^2 + c * α * (J_times_Δx ⋅ Δx)
        if armijo_condition
            break
        end

        α *= ρ
    end
    return α
end


function bt_line_search(Δx, x, composite, vars, others, xnorm; α = 1.0, ρ = 0.5, lstol=0.9, α_min = 1.0e-8)
    
    perturbed_x = @. x + α * Δx
    r       = compute_residual(composite, x, vars, others)
    rnorm   = mynorm(r, xnorm)
    
    # Iterate unless step length becomes too small
    while α > α_min
        # Apply scaled update
        perturbed_x = @. x + α * Δx
        
        # Get updated residual
        perturbed_r     = compute_residual(composite, perturbed_x, vars, others)
        perturbed_rnorm = mynorm(perturbed_r, xnorm)  

        # Check whether residual is sufficiently reduced
        if perturbed_rnorm ≤ lstol * rnorm 
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
            v += !iszero(yi) * abs(xi / yi)
        end
        return v
    end
end
