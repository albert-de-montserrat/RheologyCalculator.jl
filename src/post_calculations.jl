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

    # Walk all equations with the generic method. For elastic elements that are
    # *direct leafs* of the outer SeriesModel (IsGlobal = true equations), this
    # correctly returns xnew[self] — the effective strain-rate correction already
    # folded τ0 into the global solve, so the solved stress is the physical spring
    # stress. For elastic elements *nested inside ParallelModel branches*
    # (IsGlobal = false equations), the corresponding f-overloads now return ()
    # so those entries are skipped here and handled separately below.
    τ_direct = compute_stress_elastic(eqs, xnew, others)

    # The outer model has exactly one top-level global equation (fn =
    # compute_strain_rate, IsGlobal = true). Its solved value x[self] is the
    # shared total stress τ that all branches must carry in a series arrangement.
    # This is proven to equal the true physical stress regardless of τ0 (the
    # effective strain-rate correction cancels exactly at the global level).
    i_global = findfirst(eq -> isglobal(eq) === Val(true) && eq.fn === compute_strain_rate, eqs)
    τ_shared = xnew[eqs[i_global].self]

    # Count elastic elements that are direct leafs of the outer SeriesModel so
    # that τ0 indexing inside the branches starts at the right offset. The τ0
    # tuple is ordered: outer-series leafs first, then branch elements recursively
    # (matching global_eltype_numbering / effective_strain_rate_correction).
    n_el_leafs = count_elastic(c.leafs)

    # Reconstruct the physical spring stress for every elastic element nested
    # inside one of the ParallelModel branches, using the closed-form inversion
    # described in tensor_reduction.typ and in the docstring above.
    τ_branch = _kv_corrections_elastic_stress(c.branches, τ_shared, others, n_el_leafs)

    # Direct-leaf stresses come first, matching the τ0 tuple ordering expected
    # by the next time-step's effective_strain_rate_correction.
    return (τ_direct..., τ_branch...)
end

# -----------------------------------------------------------------------
# Branch-level elastic-stress reconstruction for SeriesModel with
# ParallelModel branches (Kelvin-Voigt / generalized Maxwell bodies)
# -----------------------------------------------------------------------
#
# Background (tensor_reduction.typ, "Mixed Kelvin-Voigt and Maxwell Body"):
#
# After the Newton solve, x[1] = τ_shared is the true physical total stress.
# For a ParallelModel branch, all sub-elements carry the same stress τ_shared,
# and the branch-local inelastic strain rate ε_p satisfies
#
#   τ_shared = 2 * η_KV * ε_p + Σ_i η_star_i * τ0_i       (*)
#
# where the sum runs over every elastic source in the branch:
#   η_star = 1               for a direct elastic leaf  (pure KV case)
#   η_star = η_eff_M / η_el  for a Maxwell sub-branch   (attenuated backstress)
#
# Inverting (*) gives ε_p, from which each spring's physical stress follows:
#   τ_spring = τ0_i + 2 * η_i * ε_p                (direct elastic leaf)
#   τ_spring = 2 * η_eff_M * (ε_p + τ0_i / (2*η_el)) (Maxwell sub-branch)
#
# All index arithmetic is resolved at *compile time* by the @generated
# functions below — the emitted code is a flat sequence of arithmetic with no
# dynamic dispatch, exactly mirroring _kv_corrections / _weighted_backstress
# in src/strain_rate_correction.jl.
# -----------------------------------------------------------------------

"""
    _kv_corrections_elastic_stress(branches, τ_shared, others, offset)

Accumulate the physical elastic spring stresses from all `ParallelModel`
branches of a `SeriesModel`.

`offset` is the number of elastic elements already accounted for by the series'
own direct leafs; `τ0[offset + k]` is the k-th elastic backstress owned by the
branches. Per-branch τ0 offsets are computed at specialisation time via
`_n_elastic_in_parallel` and baked in as literal integers.
"""
@generated function _kv_corrections_elastic_stress(branches::NTuple{Nb, Any}, τ_shared, others, offset) where {Nb}
    # Compile-time: count elastic elements in each branch so we can compute
    # the cumulative τ0 index offset for each one without any runtime bookkeeping.
    counts  = [_n_elastic_in_parallel(branches.parameters[i]) for i in 1:Nb]
    # offsets[i] = number of elastic elements consumed by branches 1 … i-1.
    offsets = cumsum([0; counts[1:(end - 1)]])

    # Build a flat statement list: one _branch_elastic_stress call per branch,
    # each with a literal offset baked in at specialisation time.
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

"""
    _branch_elastic_stress(branch::ParallelModel, τ_shared, others, el_idx_start)

Reconstruct the physical elastic spring stress for every elastic source inside
one `ParallelModel` branch.  Returns an empty tuple immediately if the branch
contains no elastic elements.

The branch-local inelastic strain rate is obtained by inverting equation (*):

    ε_p = (τ_shared - ws) / (2 * η_KV)

where `ws = Σ η_star_i * τ0_i` (scalar weighted-backstress sum) and
`η_KV` is the arithmetic-mean effective viscosity of the branch.

Each spring stress is then:
  - direct KV leaf:    `τ0_II + 2 * η * ε_p`
  - Maxwell sub-branch: `2 * η_eff_M * (ε_p + τ0_II / (2 * η_el))`

`el_idx_start` is the 0-based index of the elastic element immediately before
the first one owned by this branch (`τ0[el_idx_start + 1]` is this branch's
first backstress entry).
"""
@inline function _branch_elastic_stress(branch::ParallelModel, τ_shared, others, el_idx_start)
    # Short-circuit: no elastic content anywhere in this parallel block.
    iselastic(branch) == Val(false) && return ()

    # Build an args NamedTuple for viscosity queries. ε = 0.0 is a scalar
    # placeholder: for the linear rheologies that appear inside elastic-bearing
    # branches the result of compute_viscosity_series is independent of ε.
    # For nonlinear elements this is an approximation consistent with the level
    # of rigor already present in effective_strain_rate_correction (which also
    # evaluates viscosities at a representative ε, not the exact branch ε).
    args = merge((; ε = 0.0), others)

    # η_KV: arithmetic sum of effective viscosities of all sub-elements.
    # (See _η_KV in src/strain_rate_correction.jl for the definition.)
    η_KV = _η_KV(branch.leafs, branch.branches, args)

    # Gather per-elastic-source information (η_star, τ0_II, recovery formula
    # parameters) produced at compile time by _branch_elastic_info.
    info = _branch_elastic_info(branch.leafs, branch.branches, others.τ0, args, el_idx_start)

    # ws = Σ_i η_star_i * τ0_II_i  — the weighted backstress numerator.
    ws = sum(i -> i.η_star * i.τ0_II, info)

    # Invert equation (*) to obtain the branch-local inelastic strain rate.
    ε_p = (τ_shared - ws) / (2 * η_KV)

    # Recover the physical spring stress for each elastic source.
    return ntuple(length(info)) do k
        i = info[k]
        # Direct elastic leaf of the ParallelModel (pure KV): the spring stress
        # is simply the backstress plus the elastic loading from ε_p.
        # Maxwell SeriesModel sub-branch: use the Maxwell effective viscosity
        # and the sub-branch's elastic viscosity η_el = G * dt.
        i.direct ? i.τ0_II + 2 * i.η * ε_p : 2 * i.η_eff_M * (ε_p + i.τ0_II / (2 * i.η_el))
    end
end

"""
    _branch_elastic_info(leafs, subs, τ0, args, el_idx_start)

Compile-time–generated helper that returns a tuple of NamedTuples, one per
elastic source inside a `ParallelModel` branch, containing everything needed
by `_branch_elastic_stress` to reconstruct the physical spring stress:

| field    | meaning                                                  |
|----------|----------------------------------------------------------|
| η_star   | backstress weight: 1.0 (direct leaf) or η_eff_M/η_el    |
| τ0_II    | scalar second invariant of the previous spring stress    |
| direct   | true = direct elastic leaf, false = Maxwell sub-branch   |
| η        | elastic viscosity G*dt (direct leaf only)                |
| η_eff_M  | Maxwell effective viscosity 1/(1/η_v + 1/(G*dt)) (sub-branch) |
| η_el     | elastic viscosity G*dt of the sub-branch (sub-branch)    |

The element positions inside `leafs` and each sub-branch's leaf tuple are
found by type inspection at specialisation time (`<: AbstractElasticity`), and
the corresponding τ0 index literals are baked into the emitted code.  The
runtime body is therefore a flat sequence of multiply-adds with no branches.

`el_idx_start` carries the cumulative τ0 offset from `_kv_corrections_elastic_stress`
so that `τ0[el_idx_start + k]` is the k-th elastic backstress owned by this branch.
"""
@generated function _branch_elastic_info(leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, τ0, args, el_idx_start) where {N, Ns}
    # --- compile-time: locate elastic elements by type inspection ---
    # Positions of AbstractElasticity elements among the direct leafs of the
    # ParallelModel (these are the pure KV elastic springs).
    elastic_leaf_pos = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))

    # For each Maxwell sub-branch, positions of its elastic leaf inside *its*
    # own leaf tuple. subs.parameters[j].parameters[1] is the leaf-tuple type
    # of sub-branch j.
    sub_elastic_pos = [
        findall(Ti -> Ti <: AbstractElasticity, collect(subs.parameters[j].parameters[1].parameters))
            for j in 1:Ns
    ]

    stmts   = Any[:(info = ())]
    el_count = 0   # running count of elastic elements consumed → τ0 index

    # --- Direct elastic leafs of the ParallelModel (η_star = 1) ---
    # Each direct spring's backstress enters the weighted sum undiluted.
    for pos in elastic_leaf_pos
        el_count += 1
        push!(stmts, quote
            τ0_II = second_invariant_value(τ0[el_idx_start + $el_count])
            # η = G*dt for this elastic leaf; η_eff_M / η_el unused (direct = true).
            info = (info..., (η_star = 1.0, τ0_II = τ0_II, direct = true,
                              η = compute_viscosity(leafs[$pos], args),
                              η_eff_M = 0.0, η_el = 1.0))
        end)
    end

    # --- Maxwell SeriesModel sub-branches (η_star = η_eff_M / η_el) ---
    # The Maxwell weighting attenuates the backstress: a soft spring (small G)
    # contributes nearly all its backstress (η_star → 1); a very stiff spring
    # contributes negligibly (η_star → 0, element behaves as a pure dashpot).
    for j in 1:Ns
        for _ in sub_elastic_pos[j]
            el_count += 1
            push!(stmts, quote
                τ0_II   = second_invariant_value(τ0[el_idx_start + $el_count])
                η_eff_M = _η_eff_maxwell(subs[$j].leafs, args)   # 1/(1/η_v + 1/(G*dt))
                η_el    = _η_eff_elastic(subs[$j].leafs, args)   # G * dt
                # η = 0.0 unused for sub-branch recovery (direct = false).
                info = (info..., (η_star = η_eff_M / η_el, τ0_II = τ0_II, direct = false,
                                  η = 0.0, η_eff_M = η_eff_M, η_el = η_el))
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
