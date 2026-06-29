# Initial guess for the local solution vector

"""
    initial_guess_x(c, vars, args, others)

Return an `SVector` initial guess for the Newton solver given composite model `c`.

Internally calls [`generate_equations`](@ref) on `c` then estimates each unknown
independently.

# Arguments
- `c::AbstractCompositeModel`: composite rheology model (e.g. [`SeriesModel`](@ref), [`ParallelModel`](@ref)).
- `vars::NamedTuple`: kinematic inputs held constant during the solve
  (e.g. `(; ε = 1e-15)` for prescribed deviatoric strain rate).
- `args::NamedTuple`: initial estimates for the differentiable unknowns
  (e.g. `(; τ = 1e3)` for stress).
- `others::NamedTuple`: nondifferentiable auxiliary fields required by state
  functions, such as `dt`, `τ0`, `P0`, or grain size `d`.

# Example
```julia
c      = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
vars   = (; ε = 1e-15)
args   = (; τ = 0.0)
others = (; dt = 1e10, τ0 = (0.0,), P0 = (0.0,))
x      = initial_guess_x(c, vars, args, others)
```
"""
function initial_guess_x(c, vars, args, others)
    eqs = generate_equations(c)
    x0 = initial_guess_x(eqs, vars, args, others)
    return SA[promote(x0...)...]
end

"""
    x0 = initial_guess_x(eqs, vars, args, others)

Compute the initial guess for the local solution vector `x` given a tuple `eqs` of
`CompositeEquation`s. Each component of `x0` is estimated independently by calling
`estimate_initial_value` on the corresponding equation.

# Arguments
- `eqs::NTuple{N, CompositeEquation}`: equations generated from an `AbstractCompositeModel`,
  one per unknown in the solution vector.
- `vars::NamedTuple`: kinematic input variables held constant during the solve
  (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: current values of the differentiable unknowns (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables (e.g. `dt`, `τ0`, `P0`, `d`).

# Returns
- `x0::NTuple{N}`: tuple of scalar initial guesses, one per equation.
"""
@generated function initial_guess_x(eqs::NTuple{N, CompositeEquation}, vars, args, others) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> estimate_initial_value(eqs[i], vars, args, others)
    end
end

"""
    x_keys(c::AbstractCompositeModel)

Return the tuple of `Symbol`s identifying each entry of the solver vector `x` for
composite model `c`. The ordering matches [`generate_equations`](@ref).

Keys may repeat when multiple equations solve for the same physical unknown type
(e.g. two separate `:τ` entries in a model with nested plasticity).

# Example
```julia
c = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
x_keys(c)  # (:τ,)
```
"""
x_keys(c::AbstractCompositeModel) = x_keys(generate_equations(c))

"""
    x_keys(eqs::NTuple{N, CompositeEquation})

Return the per-equation differentiable-keyword keys for an equation tuple. Internal
overload used when the equations have already been generated.
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

"""
    estimate_initial_value(eq::CompositeEquation, vars, args, others)

Dispatch the initial-value estimation for a single equation `eq` to the appropriate
method based on the equation's kernel function (`eq.fn`):
- `compute_strain_rate` / `compute_volumetric_strain_rate` → harmonic mean over rheology components.
- `compute_stress` / `compute_pressure` → arithmetic sum over rheology components.
- any other function → returns `0`.

# Arguments
- `eq::CompositeEquation`: a single equation from the composite model, carrying the kernel
  function `eq.fn`, the tuple of rheology components `eq.rheology`, and the element
  indices `eq.el_number`.
- `vars::NamedTuple`: kinematic input variables (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: current differentiable unknowns (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables (e.g. `dt`, `τ0`, `P0`).

# Returns
- `Float64`: scalar initial-guess value for the unknown of `eq`.
"""
estimate_initial_value(eq::CompositeEquation, vars, args, others) = _estimate_initial_value(eq.fn, eq, vars, args, others)
# Fallback: unknown equation type → use 0 as the initial guess.
@inline _estimate_initial_value(::F, eq, vars, args, others) where {F} = 0
# Strain-rate-like unknowns use a harmonic-mean estimate across the element rheologies.
@inline _estimate_initial_value(::typeof(compute_volumetric_strain_rate), eq, vars, args, others) = _estimate_initial_value_harm(eq.fn, eq.rheology, eq.el_number, vars, args, others)
@inline _estimate_initial_value(::typeof(compute_strain_rate), eq, vars, args, others) = _estimate_initial_value_harm(eq.fn, eq.rheology, eq.el_number, vars, args, others)
# Stress-like unknowns use an arithmetic-sum estimate across the element rheologies.
@inline _estimate_initial_value(::typeof(compute_stress), eq, vars, args, others) = _estimate_initial_value_arith(eq.fn, eq.rheology, eq.el_number, vars, args, others)
@inline _estimate_initial_value(::typeof(compute_pressure), eq, vars, args, others) = _estimate_initial_value_arith(eq.fn, eq.rheology, eq.el_number, vars, args, others)

# Base cases for empty rheology tuples.
@inline _estimate_initial_value_harm(fn, rheology::Tuple{}, el_number, vars, args, others) = 1
@inline _estimate_initial_value_arith(fn, rheology::Tuple{}, el_number, vars, args, others) = 1

"""
    _estimate_initial_value_harm(fn, rheology, el_number, vars, args, others)

Estimate an initial value for a **strain-rate-like** unknown as the harmonic mean of
`counterpart(fn)` evaluated independently on each rheology component.

For each component `i`, history-dependent kwargs (e.g. `τ0`, `P0`) are extracted from
`others` using `el_number[i]` and merged with `args` and `vars`. The combined
`NamedTuple` is converted to scalar invariants via `tensor2invariant` before calling
the counterpart function.

# Arguments
- `fn`: the kernel state function whose counterpart is used for the estimate
  (e.g. `compute_strain_rate` → `compute_stress` is its counterpart).
- `rheology::NTuple{N, AbstractRheology}`: individual rheology components of the model.
- `el_number`: tuple of integer indices, one per component, used to retrieve
  element-local history variables from `others`.
- `vars::NamedTuple`: kinematic input variables (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: current differentiable unknowns (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables (e.g. `dt`, `τ0`, `P0`, `d`).

# Returns
- `Float64`: harmonic-mean estimate. Returns `1` when the rheology tuple is empty
  or when all component values are zero.
"""
@generated function _estimate_initial_value_harm(fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {N}
    return quote
        @inline
        sum_vals = 0.0
        Base.@nexprs $N i -> begin
            keys_hist = history_kwargs(rheology[i])
            args_local = extract_local_kwargs(others, keys_hist, el_number[i])
            args_combined = merge(args, args_local, vars)
            fn_c = counterpart(fn)
            args_invariant = tensor2invariant(args_combined)
            val = fn_c(rheology[i], args_invariant)
            # harmonic mean
            sum_vals += iszero(val) ? zero(val) : inv(val)
        end
        return iszero(sum_vals) ? one(sum_vals) : inv(sum_vals)
    end
end

"""
    _estimate_initial_value_arith(fn, rheology, el_number, vars, args, others)

Estimate an initial value for a **stress-like** unknown as the arithmetic sum of
`counterpart(fn)` evaluated independently on each rheology component.

For each component `i`, history-dependent kwargs (e.g. `τ0`, `P0`) are extracted from
`others` using `el_number[i]` and merged with `args` and `vars`. The combined
`NamedTuple` is converted to scalar invariants via `tensor2invariant` before calling
the counterpart function.

# Arguments
- `fn`: the kernel state function whose counterpart is used for the estimate
  (e.g. `compute_stress` → `compute_strain_rate` is its counterpart).
- `rheology::NTuple{N, AbstractRheology}`: individual rheology components of the model.
- `el_number`: tuple of integer indices, one per component, used to retrieve
  element-local history variables from `others`.
- `vars::NamedTuple`: kinematic input variables (e.g. `(; ε = εᵢⱼ, θ = θ)`).
- `args::NamedTuple`: current differentiable unknowns (e.g. `(; τ = τ, P = P)`).
- `others::NamedTuple`: non-differentiable auxiliary variables (e.g. `dt`, `τ0`, `P0`, `d`).

# Returns
- `Float64`: arithmetic-sum estimate. Returns `1` when the rheology tuple is empty.
"""
@generated function _estimate_initial_value_arith(fn, rheology::NTuple{N, AbstractRheology}, el_number, vars, args, others) where {N}
    return quote
        @inline
        sum_vals = 0.0
        Base.@nexprs $N i -> begin
            keys_hist = history_kwargs(rheology[i])
            args_local = extract_local_kwargs(others, keys_hist, el_number[i])
            args_combined = merge(args, args_local, vars)
            fn_c = counterpart(fn)
            args_invariant = tensor2invariant(args_combined)
            sum_vals += fn_c(rheology[i], args_invariant)
        end
        return sum_vals
    end
end

"""
    tensor2invariant(x)

Reduce tensor components to their scalar second invariant, dispatching on the type of `x`:

| Type of `x`                          | Behaviour                                                          |
|:-------------------------------------|:-------------------------------------------------------------------|
| `NTuple{N, Real}`                    | Returns `second_invariant(x...)` — scalar invariant of a flat component tuple (e.g. `(τxx, τyy, τxy)`). |
| `NTuple{N1, NTuple{N2, Real}}`       | Applies `tensor2invariant` element-wise; returns an `NTuple{N1}` of invariants, one per inner tuple. |
| `Tuple{}`                            | Returned unchanged.                                                 |
| `Number`                             | Returned unchanged; already a scalar.                              |
| `NamedTuple`                         | Applies `tensor2invariant` to every field value and reconstructs a `NamedTuple` with the same keys. |

# Arguments
- `x`: tensor data in one of the supported forms above.

# Returns
- A scalar, tuple of scalars, or `NamedTuple` of scalars depending on the input type.
"""
@inline tensor2invariant(A::Tuple{}) = A
@inline tensor2invariant(A::NTuple{N, Real}) where N = second_invariant(A...)
@generated function tensor2invariant(A::NTuple{N1, NTuple{N2, Real}}) where {N1, N2} 
    quote 
        @inline
        Base.@ntuple $N1 i -> tensor2invariant(A[i])
    end
end
@inline tensor2invariant(a::Number) = a

function tensor2invariant(x::NamedTuple)
    k = keys(x)
    v = values(x)
    invariants = ntuple(Val(length(k))) do i
        @inline 
        tensor2invariant(v[i])
    end
    (; zip(k, invariants)...)
end
