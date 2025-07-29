# Postprocessing calculations (not needed as part of local iterations)

# now update elastic stress
is_eq_elastic(::AbstractElasticity) = true
is_eq_elastic(::T) where {T} = false

"""
    τ_elastic = compute_stress_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others)

Returns stress of the elastic elements
"""
@generated function compute_stress_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others) where {N}
    return quote
        @inline
        args_all = generate_args_template(eqs, xnew, others)
        τ_elastic = Base.@ntuple $N i_eq -> begin
            args = args_all[i_eq]
            eq = eqs[i_eq]
            (; self, rheology, fn, el_number) = eq

            @inline _compute_stress_elastic(rheology, self, fn, el_number, xnew, others, args)
        end |> superflatten
        return superflatten(τ_elastic)
    end
end
compute_stress_elastic(c::AbstractCompositeModel, xnew, others) = compute_stress_elastic(generate_equations(c), xnew, others)

@generated function _compute_stress_elastic(rheology::NTuple{N, Any}, self, fn::F, el_number, xnew, others, args) where {N, F}
    return quote
        @inline
        Base.@ntuple $N i_rheo -> begin
            r = rheology[i_rheo]
            number = el_number[i_rheo]
            _compute_stress_elastic1(r, self, fn, number, xnew, others, args)
        end |> superflatten
    end
end

@inline _compute_stress_elastic1(::T, ::Vararg{Any, N}) where {T, N} = ()
function _compute_stress_elastic1(r::AbstractElasticity, self, fn::F, number, xnew, others, args) where {F}

    @inline f(::F, ::Vararg{Any, N}) where {F, N} = ()
    @inline f(::typeof(compute_strain_rate), r, others, args, number, self, xnew) = xnew[self]

    @inline function f(::typeof(compute_stress), r, others, args, number, self, xnew)
        # compute elastic stress from local strain rate
        keys_hist = history_kwargs(r)
        args_local = extract_local_kwargs(others, keys_hist, number)
        args_combined = merge(args, args_local)
        return compute_stress(r, args_combined)
    end

    return f(fn, r, others, args, number, self, xnew)
end

"""
    P_elastic = compute_pressure_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others)

Returns pressure of the elastic elements
"""
@generated function compute_pressure_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others) where {N}
    return quote
        @inline
        args_all = generate_args_template(eqs, xnew, others)
        τ_elastic = Base.@ntuple $N i_eq -> begin
            args = args_all[i_eq]
            eq = eqs[i_eq]
            (; self, rheology, fn, el_number) = eq

            @inline _compute_pressure_elastic(rheology, self, fn, el_number, xnew, others, args)
        end |> superflatten
        return superflatten(τ_elastic)
    end
end
compute_pressure_elastic(c::AbstractCompositeModel, xnew, others) = compute_pressure_elastic(generate_equations(c), xnew, others)


@generated function _compute_pressure_elastic(rheology::NTuple{N, Any}, self, fn::F, el_number, xnew, others, args) where {N, F}
    return quote
        @inline
        Base.@ntuple $N i_rheo -> begin
            r = rheology[i_rheo]
            number = el_number[i_rheo]
            _compute_pressure_elastic1(r, self, fn, number, xnew, others, args)
        end |> superflatten
    end
end

@inline _compute_pressure_elastic1(::T, ::Vararg{Any, N}) where {T, N} = ()
function _compute_pressure_elastic1(r::AbstractElasticity, self, fn::F, number, xnew, others, args) where {F}

    @inline f(::F, ::Vararg{Any, N}) where {F, N} = ()
    @inline f(::typeof(compute_volumetric_strain_rate), r, others, args, number, self, xnew) = xnew[self]

    @inline function f(::typeof(compute_pressure), r, others, args, number, self, xnew)
        # compute elastic pressure from local volumetric strain rate
        keys_hist = history_kwargs(r)
        args_local = extract_local_kwargs(others, keys_hist, number)
        args_combined = merge(args, args_local)
        return compute_pressure(r, args_combined)
    end

    return f(fn, r, others, args, number, self, xnew)
end
