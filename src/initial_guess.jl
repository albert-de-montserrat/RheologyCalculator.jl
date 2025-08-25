# Initial guess for the local solution vector

""" 
    x0 = initial_guess_x(c, vars, args, others)

Initial guess for the local solution vector `x`
"""
function initial_guess_x(c, vars, args, others)
    eqs = generate_equations(c)
    x0 = initial_guess_x(eqs, vars, args, others)
    return SVector{length(x0)}(x0)
end

@generated function initial_guess_x(eqs::NTuple{N, CompositeEquation}, vars, args, others) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> estimate_initial_value(eqs[i], vars, args, others)
    end
end

"""
    x_keys = x_keys(c::AbstractCompositeModel)
Returns the keys of the local solution vector `x` (note that they can be repeated)
"""
x_keys(c::AbstractCompositeModel) = x_keys(generate_equations(c))

"""
    x_keys = x_keys(eqs::NTuple{N,CompositeEquation})
"""
@generated function x_keys(eqs::NTuple{N, CompositeEquation}) where {N}
    return quote
        @inline
        k = Base.@ntuple $N i -> begin
            keys(differentiable_kwargs(eqs[i].fn))
        end
        superflatten(k)
    end
end

estimate_initial_value(eq::CompositeEquation, vars, args, others) = _estimate_initial_value(eq.fn, eq, vars, args, others)
@inline _estimate_initial_value(::F, eq, vars, args, others) where {F} = 0
@inline _estimate_initial_value(::typeof(compute_volumetric_strain_rate), eq, vars, args, others) = _estimate_initial_value_harm(eq.fn, eq.rheology, eq.el_number, vars, args, others)
@inline _estimate_initial_value(::typeof(compute_strain_rate), eq, vars, args, others) = _estimate_initial_value_harm(eq.fn, eq.rheology, eq.el_number, vars, args, others)
@inline _estimate_initial_value(::typeof(compute_stress), eq, vars, args, others) = _estimate_initial_value_arith(eq.fn, eq.rheology, eq.el_number, vars, args, others)
@inline _estimate_initial_value(::typeof(compute_pressure), eq, vars, args, others) = _estimate_initial_value_arith(eq.fn, eq.rheology, eq.el_number, vars, args, others)

@inline _estimate_initial_value_harm(fn, rheology::Tuple{}, el_number, vars, args, others) = 1
@inline _estimate_initial_value_arith(fn, rheology::Tuple{}, el_number, vars, args, others) = 1

@generated function _estimate_initial_value_harm(fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {N}
    return quote
        @inline
        sum_vals = 0.0
        Base.@nexprs $N i -> begin
            keys_hist = history_kwargs(rheology[i])
            args_local = extract_local_kwargs(others, keys_hist, el_number[i])
            args_combined = merge(args, args_local, vars)
            fn_c = counterpart(fn)
            val = fn_c(rheology[i], args_combined)
            # harmonic mean
            sum_vals += iszero(val) ? zero(val) : inv(val)
        end
        return iszero(sum_vals) ? one(sum_vals) : inv(sum_vals)
    end
end

@generated function _estimate_initial_value_arith(fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {N}
    return quote
        @inline
        sum_vals = 0.0
        Base.@nexprs $N i -> begin
            keys_hist = history_kwargs(rheology[i])
            args_local = extract_local_kwargs(others, keys_hist, el_number[i])
            args_combined = merge(args, args_local, vars)
            fn_c = counterpart(fn)
            sum_vals += fn_c(rheology[i], args_combined)
        end
        return sum_vals
    end
end
