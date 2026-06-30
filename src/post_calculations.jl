# Postprocessing calculations (not needed as part of local iterations)

# now update elastic stress
is_eq_elastic(::AbstractElasticity) = true
is_eq_elastic(::T) where {T} = false

"""
    compute_stress_elastic(c, xnew, others)
    compute_stress_elastic(eqs, xnew, others)

Return a flat tuple with the updated stresses of all elastic elements present in
`c` or equation tuple `eqs`.

This is a post-processing helper for advancing elastic history fields after a
successful solve. For equations written in strain-rate form, the corresponding
entry of `xnew` is returned; for equations written in stress form, the elastic
state function is evaluated from the local arguments.
"""
compute_stress_elastic(c::AbstractCompositeModel, xnew, others) = compute_stress_elastic(generate_equations(c), xnew, others)

@inline isglobal(::CompositeEquation{B}) where {B} = Val(B)

@generated function compute_stress_elastic(eqs::NTuple{N, CompositeEquation}, xnew, others) where {N}
    return quote
        @inline
        args_all = generate_args_template(eqs, xnew, others)
        τ_elastic = Base.@ntuple $N i_eq -> begin
            args = args_all[i_eq]
            eq = eqs[i_eq]
            (; self, rheology, fn, el_number) = eq
            @inline _compute_stress_elastic(rheology, self, fn, el_number, xnew, others, args, isglobal(eq))
        end |> superflatten
        return superflatten(τ_elastic)
    end
end

@generated function _compute_stress_elastic(rheology::NTuple{N, Any}, self, fn::F, el_number, xnew, others, args, eqglobal) where {N, F}
    return quote
        @inline
        Base.@ntuple $N i_rheo -> begin
            r = rheology[i_rheo]
            number = el_number[i_rheo]
            _compute_stress_elastic1(r, self, fn, number, xnew, others, args, eqglobal)
        end |> superflatten
    end
end

@inline _compute_stress_elastic1(::T, ::Vararg{Any, N}) where {T, N} = ()
function _compute_stress_elastic1(r::AbstractElasticity, self, fn::F, number, xnew, others, args, eqglobal) where {F}

    @inline f(::F, ::Vararg{Any, N}) where {F, N} = ()

    # Direct elastic leaf of the outer SeriesModel: the effective strain-rate
    # correction (see effective_strain_rate_correction / tensor_reduction.typ)
    # already absorbed τ0 into the global solve, so x[self] IS the physical
    # spring stress.
    @inline f(::typeof(compute_strain_rate), r, others, args, number, self, xnew, ::Val{true}) = xnew[self]

    # Elastic element nested inside a branch (Kelvin-Voigt / generalized
    # Maxwell sub-model): x[self] is NOT a simple ±τ0 shift of the physical
    # spring stress (see tensor_reduction.typ, mixed KV-Maxwell body — the
    # exact correction depends on the branch's full η_KV / η_star weighting,
    # not just τ0). These are reconstructed separately by
    # `compute_stress_elastic(c::SeriesModel, xnew, others)` below, using the
    # same closed-form machinery as `effective_strain_rate_correction`.
    @inline function f(::typeof(compute_strain_rate), r, others, args, number, self, xnew, ::Val{false})
        return ()
    end

    @inline function f(::typeof(compute_stress), r, others, args, number, self, xnew, ::Val{false})
        return ()
    end

    @inline function f(::typeof(compute_stress), r, others, args, number, self, xnew, ::Val{true})
        keys_hist = history_kwargs(r)
        args_local = extract_local_kwargs(others, keys_hist, number)
        args_combined = merge(args, args_local)
        τ0_II = second_invariant_value(args_local.τ0)
        return compute_stress(r, args_combined) + τ0_II
    end

    return f(fn, r, others, args, number, self, xnew, eqglobal)
end

"""
    compute_stress_elastic(c::SeriesModel, xnew::SVector, others)

Reconstruct the physical elastic spring stress for elastic elements nested
inside `ParallelModel` branches of `c` (Kelvin-Voigt and generalized Maxwell
sub-models), in addition to the direct-leaf elastic stresses already produced
by the generic `eqs`-based method.

`x[self]` for these nested equations is not a simple `±τ0` shift of the
physical spring stress; the correct relationship is the closed-form inversion
of the same η_KV / η_star weighting used by `effective_strain_rate_correction`
(tensor_reduction.typ, "Mixed Kelvin-Voigt and Maxwell Body" /
"Generalized Maxwell Body"). Given the shared top-level stress `τ_shared`
(always exactly `xnew[i]` of the model's single global equation), the local
plastic-free strain rate of a branch is

    ε_p = (τ_shared - Σ η_star_i * τ0_i) / (2 * η_KV)

from which each elastic source's physical stress is recovered as
`τ0_i + 2 * η_i * ε_p` (direct leaf) or `2 * η_eff_M_i * (ε_p + τ0_i / (2 * η_el_i))`
(Maxwell sub-branch).
"""
function compute_stress_elastic(c::SeriesModel, xnew::SVector, others)
    eqs = generate_equations(c)
    τ_direct = compute_stress_elastic(eqs, xnew, others)
    i_global = findfirst(eq -> isglobal(eq) === Val(true) && eq.fn === compute_strain_rate, eqs)
    τ_shared = xnew[eqs[i_global].self]
    n_el_leafs = count_elastic(c.leafs)
    τ_branch = _kv_corrections_elastic_stress(c.branches, τ_shared, others, n_el_leafs)
    return (τ_direct..., τ_branch...)
end

@generated function _kv_corrections_elastic_stress(branches::NTuple{Nb, Any}, τ_shared, others, offset) where {Nb}
    counts = [_n_elastic_in_parallel(branches.parameters[i]) for i in 1:Nb]
    offsets = cumsum([0; counts[1:(end - 1)]])

    stmts = Any[:(out = ())]
    for i in 1:Nb
        push!(stmts, :(out = (out..., _branch_elastic_stress(branches[$i], τ_shared, others, offset + $(offsets[i]))...)))
    end
    push!(stmts, :(out))

    return quote
        @inline
        $(stmts...)
    end
end

@inline function _branch_elastic_stress(branch::ParallelModel, τ_shared, others, el_idx_start)
    iselastic(branch) === Val(false) && return ()
    args = merge((; ε = 0.0), others)
    η_KV = _η_KV(branch.leafs, branch.branches, args)
    info = _branch_elastic_info(branch.leafs, branch.branches, others.τ0, args, el_idx_start)
    ws = sum(i -> i.η_star * i.τ0_II, info)
    ε_p = (τ_shared - ws) / (2 * η_KV)
    return ntuple(length(info)) do k
        i = info[k]
        i.direct ? i.τ0_II + 2 * i.η * ε_p : 2 * i.η_eff_M * (ε_p + i.τ0_II / (2 * i.η_el))
    end
end

@generated function _branch_elastic_info(leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, τ0, args, el_idx_start) where {N, Ns}
    elastic_leaf_pos = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))
    sub_elastic_pos = [
        findall(Ti -> Ti <: AbstractElasticity, collect(subs.parameters[j].parameters[1].parameters))
            for j in 1:Ns
    ]

    stmts = Any[:(info = ())]
    el_count = 0

    for pos in elastic_leaf_pos
        el_count += 1
        push!(stmts, quote
            τ0_II = second_invariant_value(τ0[el_idx_start + $el_count])
            info = (info..., (η_star = 1.0, τ0_II = τ0_II, direct = true, η = compute_viscosity(leafs[$pos], args), η_eff_M = 0.0, η_el = 1.0))
        end)
    end

    for j in 1:Ns
        for _ in sub_elastic_pos[j]
            el_count += 1
            push!(stmts, quote
                τ0_II = second_invariant_value(τ0[el_idx_start + $el_count])
                η_eff_M = _η_eff_maxwell(subs[$j].leafs, args)
                η_el = _η_eff_elastic(subs[$j].leafs, args)
                info = (info..., (η_star = η_eff_M / η_el, τ0_II = τ0_II, direct = false, η = 0.0, η_eff_M = η_eff_M, η_el = η_el))
            end)
        end
    end

    push!(stmts, :(info))
    return quote
        @inline
        $(stmts...)
    end
end

"""
    compute_pressure_elastic(c, xnew, others)
    compute_pressure_elastic(eqs, xnew, others)

Return a flat tuple with the updated pressures of all elastic elements present
in `c` or equation tuple `eqs`.

This is the volumetric counterpart of `compute_stress_elastic` and is intended
for updating pressure history after a successful solve.
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
