"""
    second_invariant(a)
    second_invariant(xx, yy, xy)
    second_invariant(xx, yy, zz, yz, xz, xy)

Return `a` for a scalar invariant, or the second invariant of a 2D or 3D
symmetric deviatoric tensor stored in Voigt-like component order.
"""
@inline second_invariant(a::Number) = a
@inline second_invariant(xx, yy, xy) = √((xx^2 + yy^2 + (-xx - yy)^2) / 2 + xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)
# Convenience wrappers: accept either a bare scalar or a Voigt-ordered NTuple.
@inline second_invariant_value(a::Number) = second_invariant(a)
@inline second_invariant_value(a::NTuple) = second_invariant(a...)

# -----------------------------------------------------------------------
# effective_strain_rate_correction — public entry points
# -----------------------------------------------------------------------
#
# Called once per solve step (in solve(), before the Newton loop) to convert
# the backstress histories τ0 into an additive correction on the prescribed
# strain-rate tensor, moving the elastic memory to the left-hand side:
#
#   ε_eff = ε + correction(τ0)
#
# The corrected scalar invariant εII = second_invariant(ε_eff) is then what
# the Newton solver actually sees; the elastic state functions themselves
# are τ0-free (see RheologyDefinitions.jl, tensor_reduction.typ).
# -----------------------------------------------------------------------

# Pre-solve helper for solve(): correct ONLY the direct elastic leafs of the
# outer composite using full tensor arithmetic.  ParallelModel branch corrections
# are handled implicitly inside compute_residual and must NOT be included here.
# Returns a zero correction when τ0 is absent from `others` or when the composite
# has no direct elastic leafs.
@inline function _direct_leaf_elastic_correction(c::SeriesModel, ε, others)
    hasfield(typeof(others), :τ0) || return ε .* 0
    return effective_strain_rate_correction(c.leafs, (), ε, others.τ0, others)
end
@inline _direct_leaf_elastic_correction(::AbstractCompositeModel, ε, others) = ε .* 0

"""
    effective_strain_rate_correction(c, ε, τ0, others)

Compute the effective strain-rate correction induced by previous elastic stress
history in composite `c`. Elastic elements contribute `τ0 / (2η)` using their
current effective viscosity; non-elastic elements contribute zero.
"""
effective_strain_rate_correction(c::SeriesModel, ε::NTuple, τ0::NTuple, others) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others)

# At least one elastic element exists: delegate to the (leafs, branches) decomposition.
function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε::NTuple, τ0::NTuple, others)
    effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)
end

# No elastic element anywhere in the composite: return a zero of the same shape as ε.
@inline effective_strain_rate_correction(::Val{false}, ::SeriesModel, ::NTuple{N,T}, ::NTuple{N}, ::Any) where {N,T} = zero(T)

# Scalar / non-NTuple ε overload used when ε is already a scalar invariant
# (e.g. called recursively from inside the branch correction path).
@inline effective_strain_rate_correction(c::SeriesModel, ε, τ0, others) = effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)

# Split into two independent contributions and add them:
#   1. Direct elastic leafs of the outer SeriesModel (simple Maxwell case).
#   2. Elastic elements inside ParallelModel branches (KV / generalized Maxwell).
# The τ0 tuple is ordered: leafs first (indexed 1 … n_el_leafs), then branches
# (indexed n_el_leafs+1 … end), matching global_eltype_numbering.
@inline function effective_strain_rate_correction(leafs::NTuple{N,Any}, branches::NTuple{Nb,Any}, ε, τ0::NTuple{Nτ}, others) where {N,Nb,Nτ}
    n_el_leafs = count_elastic(leafs)
    cor_leafs  = if iszero(n_el_leafs)
        ε .* 0
    else
        effective_strain_rate_correction(leafs, (), ε, τ0, others)
    end
    cor_branch = _kv_corrections(branches, ε, τ0, others, n_el_leafs)
    return cor_leafs .+ cor_branch
end

# Scan the leaf tuple for elastic elements and accumulate their corrections.
# The second argument `::Tuple{}` signals "no branches" so only leafs are handled.
# `i` is a running counter incremented only when an elastic leaf is found;
# non-elastic leafs leave i unchanged and contribute nothing (the if-guard
# prevents any τ0 access for them).  Because i is updated by a @nexprs macro
# (resolved at compile time to N literal steps), the loop has no overhead.
@generated function effective_strain_rate_correction(leafs::NTuple{N,Any}, ::Tuple{}, ε, τ0::NTuple{Nτ}, others) where {N,Nτ}
    quote
        @inline
        i = 0
        ε_elastic_cor = ε .* 0
        Base.@nexprs $N j -> begin
            i = update_correction_index(leafs[j], i)
            if i > 0
                # η = compute_viscosity(leafs[j], merge((; ε), others))
                ε_elastic_cor = ε_elastic_cor .+ effective_strain_rate_correction(leafs[j], ε, τ0[i], others, i)
            end
        end

        return ε_elastic_cor
    end
end

# Per-leaf dispatch: route through iselastic so non-elastic elements are a no-op.
@inline effective_strain_rate_correction(c::AbstractRheology, ε, τ0, others, I) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others, I)
@inline effective_strain_rate_correction(::Val{false}, c::AbstractRheology, ε, τ0, others, I) = 0

# For an elastic leaf: correction = τ0 / (2η).  This is the Maxwell backstress
# term: the strain the spring would have produced at strain rate τ0/(2G·dt)
# times dt, which must be subtracted from the total ε before solving.
# η = G·dt for an Elasticity element, so τ0/(2η) = τ0/(2G·dt).
@inline function effective_strain_rate_correction(::Val{true}, c::AbstractRheology, ε, τ0, others, I)
    η = compute_viscosity(c, merge((; ε), others))
    correction = @. τ0 / (2 * η)
    return correction
end

# Increment the elastic-element counter only when the leaf IS elastic.
# This is used instead of a runtime branch so the @nexprs unrolled loop
# can use a literal τ0 index (`τ0[i]`) without heap allocation.
@inline update_correction_index(c::AbstractRheology, I) = update_correction_index(iselastic(c), I)
@inline update_correction_index(::Val{false}, I) = I       # non-elastic: skip
@inline update_correction_index(::Val{true},  I) = I + 1  # elastic: advance τ0 cursor

# -----------------------------------------------------------------------
# Effective-viscosity aggregation helpers
# -----------------------------------------------------------------------
#
# These are used by the KV / generalized-Maxwell branch-correction machinery
# (_η_KV, _η_eff_maxwell) to combine viscosities of tuples of leaf elements.
# The harmonic mean (sum of reciprocals, then invert) is the correct formula
# for elements connected in series: the softest element dominates.
# -----------------------------------------------------------------------

# Harmonic-mean effective viscosity of a tuple of elements in *parallel*
# (stresses add → strain rates must be consistent → harmonic mean).
@generated function compute_viscosity_parallel(branches::NTuple{N,AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            η = compute_viscosity_parallel(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        η_eff = 1 / η_eff
        return η_eff
    end
end

# Harmonic-mean effective viscosity of a tuple of elements in *series*
# (strain rates add → a stiffer element relaxes faster → harmonic mean).
@generated function compute_viscosity_series(branches::NTuple{N,AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            η = compute_viscosity_series(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        η_eff = 1 / η_eff
        return η_eff
    end
end

# Public wrapper: returns a Val so callers can dispatch on the result without
# paying for a runtime branch (the Val is always resolved at compile time when
# the concrete type of `r` is known, which it always is in @generated contexts).
"""
    iselastic(r)

Return `Val(true)` when `r` is an elastic rheology or a composite containing an
elastic rheology, otherwise `Val(false)`.
"""
@inline iselastic(r::AbstractCompositeModel) = Val(_iselastic(r))
@inline iselastic(::AbstractElasticity)      = Val(true)
@inline iselastic(::AbstractRheology)        = Val(false)

# Recursive Bool-valued predicate used by the Val-returning wrappers above.
# Checks leafs first, then branches (short-circuits on first true).
@inline _iselastic(r::AbstractCompositeModel) = _iselastic(r.leafs) || _iselastic(r.branches)

# Scan a leaf tuple: return true as soon as any element is elastic.
# @nexprs unrolls the loop at compile time; early return avoids checking the rest.
@generated function _iselastic(r::NTuple{N,AbstractRheology}) where N
    quote
        @inline
        Base.@nexprs $N i -> _iselastic(r[i]) && return true
        return false
    end
end

"""
    _iselastic(r::NTuple{N, AbstractCompositeModel})

Recurse into a tuple of composite models (e.g. the `branches` field of a
`SeriesModel`) and return `true` if any of them contains an elastic element.
This extends the leaf-only `NTuple{N, AbstractRheology}` method so that
`iselastic` works correctly on nested composites.
"""
@generated function _iselastic(r::NTuple{N,AbstractCompositeModel}) where N
    quote
        @inline
        Base.@nexprs $N i -> _iselastic(r[i]) && return true
        return false
    end
end

# Tiebreaker: an empty tuple is never elastic.  Needed because Tuple{} matches
# both NTuple{N, AbstractRheology} and NTuple{N, AbstractCompositeModel} at N=0.
@inline _iselastic(::Tuple{})            = false
@inline _iselastic(::AbstractElasticity) = true
@inline _iselastic(::AbstractRheology)   = false

# -----------------------------------------------------------------------
# Generalized Maxwell / Kelvin-Voigt strain-rate correction for branches
# -----------------------------------------------------------------------
#
# Background (see tensor_reduction.typ for the full derivation):
#
# For a SeriesModel whose branches contain ParallelModel elements, each
# parallel block contributes an effective strain-rate correction beyond the
# simple Maxwell leaf correction already handled by the existing code.
#
# For one ParallelModel branch the corrected strain rate is:
#
#   ε_eff = ε + Σ_i η_star_i * τ0_i / (2 * η_KV)           (*)
#
# where the sum runs over every elastic source (leaf or sub-branch) inside
# the parallel block, and:
#
#   η_KV  = Σ η_eff_i   (sum of effective viscosities, i.e. arithmetic mean)
#   η_star = 1                                  for a direct elastic leaf
#   η_star = η_eff_M / (G * dt) = η_v / (η_v + G*dt)  for a Maxwell sub-branch
#
# η_star < 1 for Maxwell branches: a softer spring (small G) contributes more
# of its backstress; a very stiff spring makes η_star → 0 because the elastic
# strain is negligible and the element behaves like a pure dashpot.
#
# All index arithmetic and η_star expressions are resolved at *compile time*
# by the generated functions below, so the runtime code is a flat sequence of
# arithmetic with no dynamic dispatch or branching.
# -----------------------------------------------------------------------

"""
    count_elastic(r::NTuple{N, AbstractRheology})

Return, at compile time, the number of `AbstractElasticity` elements in the
leaf tuple `r`. Used to compute τ0 index offsets before processing the
`ParallelModel` branches of a `SeriesModel`.
"""
@generated function count_elastic(r::NTuple{N,AbstractRheology}) where N
    # Count by inspecting the concrete element types at specialisation time.
    n = count(T -> T <: AbstractElasticity, r.parameters)
    :($n)
end
count_elastic(::Tuple{}) = 0

"""
    _n_elastic_in_parallel(::Type{ParallelModel{L, B}})

Return, at the *type level*, the total number of elastic elements inside a
`ParallelModel`: elastic direct leafs (type `L`) plus elastic leafs of every
`SeriesModel` sub-branch in `B`.

Called only from `_kv_corrections` at specialisation time; never called at
runtime.
"""
function _n_elastic_in_parallel(::Type{ParallelModel{L,B}}) where {L,B}
    # Elastic direct leafs of this parallel block.
    n_leafs = count(T -> T <: AbstractElasticity, L.parameters)
    # For each SeriesModel sub-branch, parameters[1] is its leaf-tuple type.
    n_subs = isempty(B.parameters) ? 0 :
             sum(count(T -> T <: AbstractElasticity, S.parameters[1].parameters) for S in B.parameters)
    return n_leafs + n_subs
end

# Recursively check, at the type level, whether a nested rheology/composite
# type contains any `AbstractElasticity` element anywhere in its subtree.
_type_has_elastic(::Type{<:AbstractElasticity}) = true
_type_has_elastic(::Type{<:AbstractRheology}) = false
function _type_has_elastic(::Type{<:Union{SeriesModel{L,B},ParallelModel{L,B}}}) where {L,B}
    any(_type_has_elastic, L.parameters) || any(_type_has_elastic, B.parameters)
end

"""
    _assert_kv_nesting_supported(::Type{ParallelModel{L, B}})

The `_η_KV`/`_η_eff_maxwell`/`_weighted_backstress`/`_n_elastic_in_parallel`
formulas (see `tensor_reduction.typ`) are derived for a `ParallelModel` branch
whose `SeriesModel` sub-branches contain plain rheology leafs only -- i.e. at
most one level of Series/Parallel alternation. They do not account for a
`SeriesModel` sub-branch that itself contains a further nested `ParallelModel`.

Raise a clear error at specialisation time if such nesting contains an
elastic element, rather than silently under-counting or omitting its
backstress contribution (that element's `iselastic` still returns `true`, so
the correction would otherwise be silently incomplete instead of zero/absent).
"""
function _assert_kv_nesting_supported(::Type{ParallelModel{L,B}}) where {L,B}
    for S in B.parameters
        sub_branches = S.parameters[2]  # SeriesModel sub-branch's own `branches` field type
        for nested in sub_branches.parameters
            _type_has_elastic(nested) && error(
                "Generalized Maxwell / Kelvin-Voigt correction: an elastic element is nested " *
                "inside a ParallelModel more than one level deep inside a branch ($nested). " *
                "This is not supported by the current _η_KV / _weighted_backstress formulas, " *
                "which are only derived for a branch whose SeriesModel sub-branches contain " *
                "plain rheology leafs (see tensor_reduction.typ)."
            )
        end
    end
    return nothing
end

"""
    _kv_corrections(branches, ε, τ0, others, offset)

Accumulate the generalized Maxwell / KV effective strain-rate corrections from
all `ParallelModel` branches of a `SeriesModel`.

`offset` is the number of elastic elements already consumed by the series
leafs, so `τ0[offset + k]` is the backstress for the k-th elastic element
inside the branches.

The τ0 index for each branch is pre-computed at specialisation time (via
`_n_elastic_in_parallel`) and baked in as a literal integer, yielding
allocation-free, branch-free runtime code.
"""
@generated function _kv_corrections(
    branches::NTuple{Nb,Any}, ε, τ0, others, offset
) where Nb
    # --- compile-time: compute per-branch elastic counts and cumulative offsets ---
    foreach(_assert_kv_nesting_supported, branches.parameters)
    # How many τ0 entries does each branch own? Determined by type inspection of
    # the branch's leaf tuple (via _n_elastic_in_parallel).
    counts  = [_n_elastic_in_parallel(branches.parameters[i]) for i in 1:Nb]
    # offsets[i] = total elastic elements in branches 1 … i-1, so that
    # τ0[offset + offsets[i] + k] is the k-th entry of branch i.
    offsets = cumsum([0; counts[1:(end-1)]])

    # Emit one _kv_branch_correction call per branch, each with its τ0 start
    # index baked in as a literal integer — no runtime bookkeeping needed.
    stmts = Any[:(cor = ε .* 0)]
    for i in 1:Nb
        push!(stmts, :(cor = cor .+ _kv_branch_correction(branches[$i], ε, τ0, others, offset + $(offsets[i]))))
    end
    push!(stmts, :(cor))

    return quote
        @inline
        $(stmts...)
    end
end

"""
    _kv_branch_correction(branch::ParallelModel, ε, τ0, others, el_idx_start)

Compute the generalized Maxwell / Kelvin-Voigt effective strain-rate correction
for a single `ParallelModel` branch.  Returns zero immediately when the branch
contains no elastic elements.

The correction follows equation (*) in `tensor_reduction.typ`:

    Σ_i η_star_i * τ0_i / (2 * η_KV)

`el_idx_start` is the 1-based index of the elastic element immediately before
the first elastic element owned by this branch (i.e. `τ0[el_idx_start + 1]`
is this branch's first backstress entry).
"""
@inline function _kv_branch_correction(branch::ParallelModel, ε, τ0, others, el_idx_start)
    # Short-circuit: no elastic elements anywhere in this parallel block.
    iselastic(branch) == Val(false) && return ε .* 0
    # Reduce ε to its scalar second invariant for viscosity queries.
    # All viscosity functions operate on the invariant, not the full tensor.
    εII  = second_invariant_value(ε)
    args = merge((; ε = εII), others)
    # η_KV: arithmetic sum of effective viscosities of all sub-elements
    # (viscous leafs + Maxwell sub-branches).  This is the denominator of (*).
    η_KV = _η_KV(branch.leafs, branch.branches, args)
    # ws = Σ_i η_star_i * τ0_i — the weighted backstress numerator of (*).
    # Same shape as ε (tensor), because τ0 entries are also stored as tensors.
    ws = _weighted_backstress(branch.leafs, branch.branches, ε, τ0, args, el_idx_start)
    # Final correction tensor: ws / (2 * η_KV) broadcasted element-wise.
    return @. ws / (2 * η_KV)
end

"""
    _η_KV(leafs, subs, args)

Compute the Kelvin-Voigt effective viscosity of a parallel block:

    η_KV = Σ_j η_j  (viscous leafs)  +  Σ_k η_eff_M_k  (Maxwell sub-branches)

where `η_eff_M = 1 / (1/η_v + 1/(G*dt))` is the harmonic-mean effective
viscosity of each Maxwell `SeriesModel` sub-branch.
"""
@generated function _η_KV(leafs::NTuple{N,AbstractRheology}, subs::NTuple{Ns,Any}, args) where {N,Ns}
    quote
        @inline
        η = 0.0
        # Arithmetic sum over direct leafs of the ParallelModel.
        # For a linear viscous leaf this is just η_leaf; for an elastic leaf
        # (G*dt) it also contributes to the KV stiffness.
        Base.@nexprs $N  i -> η += compute_viscosity_series(leafs[i], args)
        # Each Maxwell SeriesModel sub-branch contributes its harmonic-mean
        # effective viscosity η_eff_M = (η_v * G*dt)/(η_v + G*dt).
        Base.@nexprs $Ns j -> η += _η_eff_maxwell(subs[j].leafs, args)
        η
    end
end

# Specialisation for a pure KV block with no Maxwell sub-branches.
@generated function _η_KV(leafs::NTuple{N,AbstractRheology}, ::Tuple{}, args) where N
    quote
        @inline
        η = 0.0
        Base.@nexprs $N i -> η += compute_viscosity_series(leafs[i], args)
        η
    end
end

"""
    _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args)

Effective viscosity of a Maxwell `SeriesModel` branch: the harmonic mean of
the viscosities of its constituent leaf elements.

    η_eff_M = 1 / Σ_i (1 / η_i)

For a two-element branch `(viscous, elastic)` this reduces to the standard
Maxwell formula `η_v * G * dt / (η_v + G * dt)`.
"""
@generated function _η_eff_maxwell(leafs::NTuple{N,AbstractRheology}, args) where N
    quote
        @inline
        inv_η = 0.0
        # Accumulate 1/η for each leaf; the harmonic mean (inverse of sum of
        # inverses) is the correct effective viscosity for elements in series.
        Base.@nexprs $N i -> inv_η += inv(compute_viscosity_series(leafs[i], args))
        inv(inv_η)
    end
end

"""
    _η_eff_elastic(leafs, args)

Return the effective viscosity of the elastic element inside a Maxwell
sub-branch leaf tuple (= `G * dt`).  The elastic element is identified at
*compile time* by type inspection, so the generated code is a single
`compute_viscosity` call with a literal index.

Returns `0.0` if no elastic element is found (should not happen for a Maxwell
branch, but is safe for non-elastic tuples).
"""
@generated function _η_eff_elastic(leafs::T, args) where T
    # Resolve the elastic element's position at specialisation time by scanning
    # the concrete leaf types. The emitted code is a single compute_viscosity
    # call with a literal index — no runtime search.
    idx = findfirst(Ti -> Ti <: AbstractElasticity, collect(T.parameters))
    idx === nothing && return :(0.0)  # safe fallback; should not occur for a Maxwell branch
    :(compute_viscosity(leafs[$idx], args))
end

"""
    _weighted_backstress(leafs, subs, ε, τ0, args, el_idx_start)

Compute the weighted backstress numerator `Σ_i η_star_i * τ0_i` for a
`ParallelModel` branch.

Two classes of elastic sources contribute:
- **Direct elastic leafs** of the `ParallelModel`: `η_star = 1`.  These
  correspond to the simple Kelvin-Voigt case where the elastic element is a
  direct parallel element (backstress enters undiluted).
- **Maxwell `SeriesModel` sub-branches**: `η_star = η_eff_M / η_el`, where
  `η_eff_M` is the harmonic-mean effective viscosity of the sub-branch and
  `η_el = G * dt`.  This is the generalized Maxwell weighting: a softer spring
  (small G) makes `η_star → 1`; a very stiff spring makes `η_star → 0`.

All τ0 index literals and η_star computations are resolved at *compile time*
(via `@generated`), so the emitted code is a flat sequence of multiply-adds.
`el_idx_start` carries the τ0 offset inherited from the outer `SeriesModel`
leaf count.
"""
@generated function _weighted_backstress(
    leafs::NTuple{N,AbstractRheology}, subs::NTuple{Ns,Any}, ε, τ0, args, el_idx_start
) where {N,Ns}
    # --- compile-time: locate elastic elements by type inspection ---
    # Positions of elastic elements among the direct leafs of the ParallelModel.
    elastic_leaf_pos = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))
    # For each sub-branch, positions of elastic leafs in *its* leaf tuple.
    # subs.parameters[j].parameters[1] is the leaf-tuple type of sub-branch j.
    sub_elastic_pos = [
        findall(Ti -> Ti <: AbstractElasticity, collect(subs.parameters[j].parameters[1].parameters))
        for j in 1:Ns
    ]

    stmts    = Any[:(ws = ε .* 0)]  # tensor accumulator, same shape as ε
    el_count = 0                     # compile-time cursor into the τ0 tuple

    # --- Direct elastic leafs of the ParallelModel: η_star = 1 ---
    # The backstress enters the weighted sum undiluted (pure KV case).
    # Each iteration bakes a literal τ0 index into the emitted code.
    for _ in elastic_leaf_pos
        el_count += 1
        push!(stmts, :(ws = ws .+ τ0[el_idx_start + $el_count]))
    end

    # --- Maxwell SeriesModel sub-branches: η_star = η_eff_M / η_el ---
    # Each Maxwell sub-branch attenuates its backstress by the ratio of its
    # harmonic-mean effective viscosity to the elastic viscosity (G * dt).
    for j in 1:Ns
        for _ in sub_elastic_pos[j]
            el_count += 1
            push!(stmts, quote
                η_eff_M = _η_eff_maxwell(subs[$j].leafs, args)   # 1/(1/η_v + 1/(G*dt))
                η_el = _η_eff_elastic(subs[$j].leafs, args)   # G * dt
                η_star = η_eff_M / η_el                         # η_v / (η_v + G*dt) ∈ (0,1)
                ws = ws .+ η_star .* τ0[el_idx_start+$el_count]
            end)
        end
    end

    push!(stmts, :(ws))
    return quote
        @inline
        $(stmts...)
    end
end

# -----------------------------------------------------------------------
# Implicit elastic correction — added to the Newton residual (Option 3)
# -----------------------------------------------------------------------
#
# Rather than pre-correcting the input strain rate before the Newton loop,
# these functions subtract the elastic backstress correction directly from
# the global-equation residual at each Newton iteration:
#
#   R[1]  -=  correction(x)
#
# where  correction(x) = Σ direct-leaf τ0/(2η)
#                       + Σ_branches ws(x[branch_eq]) / (2 η_KV(x[branch_eq]))
#
# Because the branch terms depend on x (through x[branch_eq_idx], the branch
# strain rate), ForwardDiff automatically differentiates through them, giving
# an exact Jacobian for any viscosity law.
#
# For LINEAR viscosities η = const and results are numerically identical to
# the old pre-correction approach.  For nonlinear viscosities in branches
# the Newton iterations now converge to the true fully-implicit solution.
# -----------------------------------------------------------------------

"""
    _implicit_elastic_correction(c::SeriesModel, eqs, x, others)

Scalar elastic backstress correction to subtract from the global residual.
Combines the constant direct-leaf contribution and the x-dependent branch
contribution (which is differentiated through by ForwardDiff).
"""
function _implicit_elastic_correction(c::SeriesModel, eqs, x::SVector{N, T}, others) where {N, T}
    # τ0 may be absent when compute_residual is called directly without elastic
    # history (e.g. test code or initialisation).  No correction in that case.
    hasfield(typeof(others), :τ0) || return zero(T)
    τ0         = others.τ0
    n_el_leafs = count_elastic(c.leafs)
    # Direct leaf (Maxwell) corrections are pre-applied in solve() with full tensor
    # arithmetic.  Only the x-dependent branch corrections belong here.
    return _kv_implicit_corrections_scalar(c.branches, eqs, x, τ0, others, n_el_leafs)
end

# Constant correction from elastic elements that are direct leafs of the outer
# SeriesModel (simple Maxwell: each contributes τ0_II / (2G·dt)).
# Elastic leaf positions and τ0 indices are resolved at specialisation time.
@generated function _direct_leaf_correction_scalar(leafs::NTuple{N, Any}, τ0, others) where {N}
    elastic_positions = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))
    stmts = Any[:(cor = 0.0)]
    for (count, pos) in enumerate(elastic_positions)
        push!(stmts, quote
            η    = compute_viscosity(leafs[$pos], others)
            cor += second_invariant_value(τ0[$count]) / (2 * η)
        end)
    end
    push!(stmts, :(cor))
    return quote @inline; $(stmts...) end
end

# Return the position within `positions` (in order) of the i-th entry whose
# `mask` is true. `positions`/`mask` are small homogeneous tuples (Int/Bool),
# so this is fully type-stable regardless of how heterogeneous the underlying
# `eqs::NTuple{N, CompositeEquation}` tuple is.
@inline function _nth_true_position(mask::NTuple{M, Bool}, positions::NTuple{M, Int}, i::Int) where {M}
    c = 0
    for j in 1:M
        if mask[j]
            c += 1
            c == i && return positions[j]
        end
    end
    error("_nth_true_position: found fewer than $i matching branch equations (found $c)")
end

# Accumulate the implicit KV/Maxwell branch corrections across all branches.
#
# A top-level branch's *own* compute_stress equation is uniquely identified by
# `eq.parent == eqs[1].self` (the outer SeriesModel's own equation is always
# the first equation emitted, self = 1) combined with `fn === compute_stress`:
# only a branch's own equation is a *direct* child of the outer equation;
# anything nested further inside that branch (e.g. a ParallelModel nested more
# than one level deep, or another branch's own sub-structure) has `.parent`
# pointing at its own enclosing equation instead. This correctly disambiguates
# branches from each other and from their own nested sub-structure regardless
# of nesting depth or whether rheology types collide anywhere in the tree —
# unlike matching on rheology type alone, which can be fooled by a *different*
# branch's nested equation sharing the same leaf type.
#
# `.parent`/`.self` are runtime Int64 fields (not part of `CompositeEquation`'s
# type), so this match can't be resolved at `@generated` specialisation time;
# only the `fn === compute_stress` pre-filter can (CompositeEquation{IsGlobal,
# T, F, R, RT}: F is parameter index 3). The pre-filtered candidate positions
# are baked in as literals, and the final `.parent` comparisons run once per
# branch at call time via `_nth_true_position` above.
@generated function _kv_implicit_corrections_scalar(
    branches::NTuple{Nb, Any}, eqs::NTuple{N, Any}, x, τ0, others, offset
) where {Nb, N}
    foreach(_assert_kv_nesting_supported, branches.parameters)
    stress_positions = Tuple(k for k in 1:N if eqs.parameters[k].parameters[3] === typeof(compute_stress))
    mask_expr = Expr(:tuple, (:(eqs[$k].parent == outer_self) for k in stress_positions)...)

    counts  = [_n_elastic_in_parallel(branches.parameters[i]) for i in 1:Nb]
    offsets = cumsum([0; counts[1:(end - 1)]])

    stmts = Any[:(cor = 0.0), :(outer_self = eqs[1].self), :(mask = $mask_expr)]
    for i in 1:Nb
        push!(stmts, quote
            bpos = _nth_true_position(mask, $stress_positions, $i)
            cor += _kv_implicit_branch_correction_scalar(
                branches[$i], x[bpos], τ0, others, offset + $(offsets[i])
            )
        end)
    end
    push!(stmts, :(cor))
    return quote @inline; $(stmts...) end
end

_kv_implicit_corrections_scalar(::Tuple{}, ::NTuple{N, Any}, ::Any, ::Any, ::Any, ::Any) where {N} = 0.0

# Per-branch implicit correction: ws(ε_branch) / (2 η_KV(ε_branch)).
# ε_branch = x[branch_eq_idx] is the branch's local strain rate from the
# current Newton iterate — correct for both linear and nonlinear viscosities.
@inline function _kv_implicit_branch_correction_scalar(branch::ParallelModel, ε_branch, τ0, others, el_idx_start)
    iselastic(branch) == Val(false) && return zero(ε_branch)
    args = merge((; ε = ε_branch), others)
    η_KV = _η_KV(branch.leafs, branch.branches, args)
    ws   = _weighted_backstress_scalar(branch.leafs, branch.branches, τ0, args, el_idx_start)
    return ws / (2 * η_KV)
end

# Scalar version of _weighted_backstress: returns Σ_i η_star_i * τ0_II_i.
# Structurally identical to _weighted_backstress but accumulates a scalar rather
# than a tensor; τ0 indices and η_star expressions are compiled-in as literals.
@generated function _weighted_backstress_scalar(
    leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, τ0, args, el_idx_start
) where {N, Ns}
    elastic_leaf_pos = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))
    sub_elastic_pos  = [
        findall(Ti -> Ti <: AbstractElasticity, collect(subs.parameters[j].parameters[1].parameters))
        for j in 1:Ns
    ]

    stmts    = Any[:(ws = 0.0)]
    el_count = 0

    for _ in elastic_leaf_pos
        el_count += 1
        push!(stmts, :(ws += second_invariant_value(τ0[el_idx_start + $el_count])))
    end

    for j in 1:Ns
        for _ in sub_elastic_pos[j]
            el_count += 1
            push!(stmts, quote
                η_eff_M = _η_eff_maxwell(subs[$j].leafs, args)
                η_el    = _η_eff_elastic(subs[$j].leafs, args)
                ws     += (η_eff_M / η_el) * second_invariant_value(τ0[el_idx_start + $el_count])
            end)
        end
    end

    push!(stmts, :(ws))
    return quote @inline; $(stmts...) end
end
