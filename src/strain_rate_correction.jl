"""
    second_invariant(a)
    second_invariant(xx, yy, xy)
    second_invariant(xx, yy, zz, yz, xz, xy)

Return `a` for a scalar invariant, or the second invariant of a 2D or 3D
symmetric deviatoric tensor stored in Voigt-like component order.
"""
@inline second_invariant(a::Number) = a
@inline second_invariant(xx, yy, xy) = √((xx^2 + yy^2 + (-xx - yy)^2) / 2 +  xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)
@inline second_invariant_value(a::Number) = second_invariant(a)
@inline second_invariant_value(a::NTuple) = second_invariant(a...)

"""
    effective_strain_rate_correction(c, ε, τ0, others)

Compute the effective strain-rate correction induced by previous elastic stress
history in composite `c`. Elastic elements contribute `τ0 / (2η)` using their
current effective viscosity; non-elastic elements contribute zero.
"""
effective_strain_rate_correction(c::SeriesModel, ε::NTuple, τ0::NTuple, others) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others)

# @generated function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε::NTuple, τ0::NTuple{N}, others) where N
#     quote
#         @inline 
#         Base.@ntuple $N i -> effective_strain_rate_correction(c, ε, τ0[i], others)
#     end
# end
function effective_strain_rate_correction(::Val{true}, c::SeriesModel, ε::NTuple, τ0::NTuple, others)
    effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)
end

@inline effective_strain_rate_correction(::Val{false}, c::SeriesModel, ε::NTuple{N, T}, τ0::NTuple{N}, others) where {N,T} = zero(T)

@inline effective_strain_rate_correction(c::SeriesModel, ε, τ0, others) = effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)

@inline function effective_strain_rate_correction(leafs::NTuple{N, Any}, branches::NTuple{Nb, Any}, ε, τ0::NTuple{Nτ}, others) where {N, Nb, Nτ}
    cor_leafs  = effective_strain_rate_correction(leafs, (), ε, τ0, others)
    n_el_leafs = count_elastic(leafs)
    cor_branch = _kv_corrections(branches, ε, τ0, others, n_el_leafs)
    return cor_leafs .+ cor_branch
end

@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0::NTuple{Nτ}, others) where {N, Nτ}
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

@inline effective_strain_rate_correction(c::AbstractRheology, ε, τ0, others, I) = effective_strain_rate_correction(iselastic(c), c, ε, τ0, others, I)
@inline effective_strain_rate_correction(::Val{false}, c::AbstractRheology, ε, τ0, others, I) = 0

@inline function effective_strain_rate_correction(::Val{true}, c::AbstractRheology, ε, τ0, others, I)
    η = compute_viscosity(c, merge((; ε), others))
    correction = @. τ0 / (2 * η)
    return correction
end

@inline update_correction_index(c::AbstractRheology, I) = update_correction_index(iselastic(c), I)
@inline update_correction_index(::Val{false}, I) = I
@inline update_correction_index(::Val{true}, I)  = I + 1

@generated function compute_viscosity_parallel(branches::NTuple{N, AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            # compute local effect on strainrate tensor due to old elastic stresses
            η = compute_viscosity_parallel(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        η_eff = 1 / η_eff
        return η_eff
    end
end

@generated function compute_viscosity_series(branches::NTuple{N, AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            # compute local effect on strainrate tensor due to old elastic stresses
            η = compute_viscosity_series(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        η_eff = 1 / η_eff
        return η_eff
    end
end

"""
    iselastic(r)

Return `Val(true)` when `r` is an elastic rheology or a composite containing an
elastic rheology, otherwise `Val(false)`.
"""
@inline iselastic(r::AbstractCompositeModel) = Val(_iselastic(r))
@inline iselastic(::AbstractElasticity)      = Val(true)
@inline iselastic(::AbstractRheology)        = Val(false)

@inline _iselastic(r::AbstractCompositeModel) = _iselastic(r.leafs) || _iselastic(r.branches)

@generated function _iselastic(r::NTuple{N, AbstractRheology}) where N
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
@generated function _iselastic(r::NTuple{N, AbstractCompositeModel}) where N
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
@generated function count_elastic(r::NTuple{N, AbstractRheology}) where N
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
function _n_elastic_in_parallel(::Type{ParallelModel{L, B}}) where {L, B}
    # Elastic direct leafs of this parallel block.
    n_leafs = count(T -> T <: AbstractElasticity, L.parameters)
    # For each SeriesModel sub-branch, parameters[1] is its leaf-tuple type.
    n_subs  = isempty(B.parameters) ? 0 :
        sum(count(T -> T <: AbstractElasticity, S.parameters[1].parameters) for S in B.parameters)
    return n_leafs + n_subs
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
    branches::NTuple{Nb, Any}, ε, τ0, others, offset
) where Nb
    # --- compile-time: compute per-branch elastic counts and cumulative offsets ---
    counts  = [_n_elastic_in_parallel(branches.parameters[i]) for i in 1:Nb]
    # offsets[i] = number of elastic elements consumed by branches 1 … i-1
    offsets = cumsum([0; counts[1:end-1]])

    # Build the flat statement list: one call per branch with a literal offset.
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
    # Reduce ε to its second invariant for viscosity queries (scalar path).
    εII  = second_invariant_value(ε)
    args = merge((; ε = εII), others)
    # η_KV: arithmetic sum of all effective viscosities in this parallel block.
    η_KV = _η_KV(branch.leafs, branch.branches, args)
    # Weighted backstress numerator: Σ η_star_i * τ0_i.
    ws   = _weighted_backstress(branch.leafs, branch.branches, ε, τ0, args, el_idx_start)
    return @. ws / (2 * η_KV)
end

"""
    _η_KV(leafs, subs, args)

Compute the Kelvin-Voigt effective viscosity of a parallel block:

    η_KV = Σ_j η_j  (viscous leafs)  +  Σ_k η_eff_M_k  (Maxwell sub-branches)

where `η_eff_M = 1 / (1/η_v + 1/(G*dt))` is the harmonic-mean effective
viscosity of each Maxwell `SeriesModel` sub-branch.
"""
@generated function _η_KV(leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, args) where {N, Ns}
    quote
        @inline
        η = 0.0
        Base.@nexprs $N  i -> η += compute_viscosity_series(leafs[i], args)
        Base.@nexprs $Ns j -> η += _η_eff_maxwell(subs[j].leafs, args)
        η
    end
end

# Specialisation for an empty sub-branch tuple (no Maxwell branches, only leafs).
@generated function _η_KV(leafs::NTuple{N, AbstractRheology}, ::Tuple{}, args) where N
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
@generated function _η_eff_maxwell(leafs::NTuple{N, AbstractRheology}, args) where N
    quote
        @inline
        inv_η = 0.0
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
    # Resolve the elastic element's position at specialisation time.
    idx = findfirst(Ti -> Ti <: AbstractElasticity, collect(T.parameters))
    idx === nothing && return :(0.0)
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
    leafs::NTuple{N, AbstractRheology}, subs::NTuple{Ns, Any}, ε, τ0, args, el_idx_start
) where {N, Ns}
    # --- compile-time: locate elastic elements by type inspection ---
    # Positions of elastic elements among the direct leafs of the ParallelModel.
    elastic_leaf_pos = findall(Ti -> Ti <: AbstractElasticity, collect(leafs.parameters))
    # For each sub-branch, positions of elastic leafs in *its* leaf tuple.
    # subs.parameters[j].parameters[1] is the leaf-tuple type of sub-branch j.
    sub_elastic_pos  = [
        findall(Ti -> Ti <: AbstractElasticity, collect(subs.parameters[j].parameters[1].parameters))
        for j in 1:Ns
    ]

    stmts    = Any[:(ws = ε .* 0)]   # accumulator, same shape as ε
    el_count = 0                      # running count of elastic elements consumed

    # --- elastic leafs of the ParallelModel: η_star = 1 ---
    # The backstress enters without attenuation; this is the pure KV case.
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
                η_el    = _η_eff_elastic(subs[$j].leafs, args)   # G * dt
                η_star  = η_eff_M / η_el                         # η_v / (η_v + G*dt) ∈ (0,1)
                ws      = ws .+ η_star .* τ0[el_idx_start + $el_count]
            end)
        end
    end

    push!(stmts, :(ws))
    return quote
        @inline
        $(stmts...)
    end
end
