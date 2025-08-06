function compute_residual(c, x::SVector{N, T}, vars, others) where {N, T}

    eqs = generate_equations(c)
    @assert length(eqs) == length(x)
    args_all = generate_args_template(eqs, x, others)

    # evaluates the self-components of the residual
    residual1 = evaluate_state_functions(eqs, args_all, others)
    residual2 = add_children(residual1, x, eqs)
    residual3 = subtract_parent(residual2, x, eqs, vars)

    return SA[residual3...]
end

function compute_residual(c, x::SVector{N, T}, vars, others, ::Int64, ::Int64) where {N, T}

    eqs = generate_equations(c)
    args_all = first(generate_args_template(eqs, x, others))

    # evaluates the self-components of the residual
    eq = first(eqs)
    residual1 = evaluate_state_function(eq, args_all, others)
    residual2 = add_children(residual1, x, eq)
    residual3 = subtract_parent(residual2, x, eq, vars)

    return residual3
end

@inline function evaluate_state_function(eq::CompositeEquation, args, others)
    (; fn, rheology, el_number) = eq
    return evaluate_state_function(fn, rheology, args, others, el_number)
end

@inline evaluate_state_function(fn::F, rheology::Tuple{}, args, others, el_number) where {F} = 0

@generated function evaluate_state_function(fn::F, rheology::NTuple{N, AbstractRheology}, args, others, el_number) where {N, F}
    return quote
        @inline
        vals = Base.@ntuple $N i -> begin
            keys_hist = history_kwargs(rheology[i])
            args_local = extract_local_kwargs(others, keys_hist, el_number[i])
            args_combined = merge(args, args_local)
            fn(rheology[i], args_combined)
        end
        sum(vals)
    end
end

evaluate_state_function(fn::F, rheology::Tuple{}, args, others) where {F} = 0.0e0

# @inline evaluate_state_functions(eqs::NTuple{N, CompositeEquation}, args) where N = promote(ntuple(i -> evaluate_state_function(eqs[i], args[i]), Val(N))...)
@generated function evaluate_state_functions(eqs::NTuple{N, CompositeEquation}, args, others) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> evaluate_state_function(eqs[i], args[i], others)
    end
end

add_child(::SVector, ::Tuple{}) = 0.0e0

@generated function add_child(x::SVector{M}, child::NTuple{N}) where {M, N}
    return quote
        @inline
        v = Base.@ntuple $N i -> begin
            x[child[i]]
        end
        sum(v)
    end
end

@generated function add_children(residual::NTuple{N, Any}, x::SVector{N}, eqs::NTuple{N, CompositeEquation}) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> residual[i] + add_child(x, eqs[i].child)
    end
end

function add_children(residual::Number, x::SVector, eq::CompositeEquation)
    return residual + add_child(x, eq.child)
end

# if global, subtract the variables
@inline subtract_parent(::SVector, eq::CompositeEquation{true}, vars) = vars[eq.ind_input]
@inline subtract_parent(x::SVector, eq::CompositeEquation{false}, ::NamedTuple) = x[eq.parent]
# exception for lambda
@inline subtract_parent(x::SVector, eq::CompositeEquation{false, T, typeof(compute_lambda)}, ::NamedTuple) where {T} = 0 # x[eq.self]

@generated function subtract_parent(residual::NTuple{N, Any}, x, eqs::NTuple{N, CompositeEquation}, vars) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> residual[i] - subtract_parent(x, eqs[i], vars)
    end
end

subtract_parent(residual::Number, x::SVector, eq::CompositeEquation, vars) = residual - subtract_parent(x, eq, vars)


