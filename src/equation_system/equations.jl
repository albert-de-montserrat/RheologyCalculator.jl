"""
    CompositeEquation

Internal representation of one residual equation generated from a composite
rheology model.

Fields store the equation index (`self`), parent equation index, child equation
indices, state function `fn`, rheology tuple, input-variable index, and
type-local element numbering for history lookup.
"""
struct CompositeEquation{IsGlobal, T, F, R, RT}
    parent::Int64       # i-th element of x to be substracted
    child::T            # i-th element of x to be added
    self::Int64         # equation number
    fn::F               # state function
    rheology::R
    ind_input::Int64
    el_number::RT

    function CompositeEquation(parent::Int64, child::T, self::Int64, fn::F, rheology::R, ind_input, ::Val{B}, el_number::RT) where {T, F, R, B, RT}
        @assert B isa Bool
        return new{B, T, F, R, RT}(parent, child, self, fn, rheology, ind_input, el_number)

    end
end

"""
    generate_equations(c::AbstractCompositeModel)

Generate the tuple of `CompositeEquation`s that defines the residual system for
composite model `c`. The length of this tuple is the expected length of the
solver vector `x`.
"""
@inline generate_equations(c::AbstractCompositeModel) = generate_equations(c, global_series_functions(c))

@generated function generate_equations(c::AbstractCompositeModel, fns::NTuple{N, Any}) where {N}
    return quote
        @inline
        iparent = 0
        iself = 0
        isGlobal = Val(true)
        # global element numbering (to be followed)
        el_num = global_eltype_numbering(c)
        eqs = Base.@ntuple $N i -> begin
            eqs = generate_equations(c, fns[i], i, isGlobal, isvolumetric(c), el_num; iparent = iparent, iself = iself)
            iself = eqs[end].self
            iparent = 0
            eqs
        end
        superflatten(eqs)
    end
end

function generate_equations(c::AbstractCompositeModel, fns_own_global::F, ind_input, ::Val{B}, ::Val, el_num; iparent::Int64 = 0, iself::Int64 = 0) where {F, B}

    @inline foo(::NTuple{N, Any}) where {N} = Val(N)

    iself_ref = Ref{Int64}(iself)
    (; branches, leafs) = c
    local_el = el_num[1]

    _, fns_own_local = get_own_functions(c)
    # fns_branches_global,_ = get_own_functions(branches)

    nown = 1 # length(fns_own_global)
    Nlocal = foo(fns_own_local)

    ilocal_childs = generate_ilocal_childs(iself, nown, fns_own_local)

    # The equations contributed by each branch are not known yet: a branch's
    # equation count depends on its own internal structure (e.g. a nested
    # SeriesModel contributes more than one equation), not merely on the
    # number of branches. The branch children are therefore left empty here
    # and patched in below, once the branches' own equations have actually
    # been generated (see attach_parallel_children).
    isGlobal = Val(B)
    global_eqs0 = add_global_equations(iparent, ilocal_childs, iself_ref, fns_own_global, leafs, ind_input, isGlobal, Val(1), local_el)
    local_eqs = add_local_equations(global_eqs0.self, (), iself_ref, fns_own_local, fns_own_global, leafs, Nlocal, local_el)
    # need to correct the children of the global equations in absence of local equations (i.e. remove them)
    global_eqs1 = correct_children(global_eqs0, local_eqs)

    fn = counterpart(fns_own_global)
    parallel_eqs = generate_equations_unroller(branches, fn, el_num, global_eqs1, iself_ref)
    global_eqs = attach_parallel_children(global_eqs1, parallel_eqs)

    # return  global_eqs, parallel_eqs
    return (global_eqs, local_eqs..., parallel_eqs...) |> superflatten
end

@generated function generate_equations_unroller(branches::NTuple{N, Any}, fn::F, el_num, global_eqs, iself_ref) where {N, F}
    return quote
        @inline
        Base.@ntuple $N i -> begin
            b = branches[i]
            num = el_num[2][i]
            isvol = isvolumetric(b)
            eqs = generate_equations(b, fn, 0, Val(false), isvol, num; iparent = global_eqs.self, iself = iself_ref[])
            # generate_equations creates its own *local* Ref for `b`'s subtree, so
            # the running counter must be threaded back here explicitly, otherwise
            # every sibling branch starts numbering its own equations from the same
            # base and their `.self` values collide (see CLAUDE.md: equation
            # position is assumed to equal `.self` throughout the residual
            # assembly). A branch can legitimately contribute zero equations for
            # a given global function (e.g. a non-volumetric branch during the
            # volumetric pass), in which case the counter is left untouched.
            isempty(eqs) || (iself_ref[] = eqs[end].self)
            eqs
        end
    end
end

generate_equations_unroller(::Tuple{}, ::F, el_num, global_eqs, iself_ref) where {F} = ()

@generated function generate_ilocal_childs(iself, nown, ::NTuple{N, Any}) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> iself + nown + i
    end
end

# A branch's own equations (`parallel_eqs[i]`) are generated after the parent
# equation, so the parent's children referencing them cannot be computed
# up-front from branch counts/offsets alone (a branch may contribute more
# than one equation, e.g. a nested SeriesModel). Instead, once the branches'
# equations exist, we patch the `.self` of each branch's top-level equation
# into the parent's `child` tuple. A branch that contributed no equations
# (e.g. not applicable to the current global function) is simply omitted.
@inline first_self(::Tuple{}) = ()
@inline first_self(eqs::Tuple) = (first(eqs).self,)

@generated function collect_parallel_children(parallel_eqs::NTuple{N, Any}) where {N}
    return quote
        @inline
        tops = Base.@ntuple $N i -> first_self(parallel_eqs[i])
        superflatten(tops)
    end
end

function attach_parallel_children(eqs::CompositeEquation{IsGlobal, T, F, R, RT}, parallel_eqs) where {IsGlobal, T, F, R, RT}
    (; parent, child, self, fn, rheology, ind_input, el_number) = eqs
    children = (child..., collect_parallel_children(parallel_eqs)...)
    return CompositeEquation(parent, children, self, fn, rheology, ind_input, Val(IsGlobal), el_number)
end


correct_children(eqs::CompositeEquation, ::NTuple{N, CompositeEquation}) where {N} = eqs

function correct_children(eqs::CompositeEquation{IsGlobal, T, F, R, RT}, ::Tuple{Tuple{}}) where {IsGlobal, T, F, R, RT}
    (; parent, self, fn, rheology, ind_input, el_number) = eqs
    return CompositeEquation(parent, (), self, fn, rheology, ind_input, Val(IsGlobal), el_number)
end

# eliminate equations
for fn in (:compute_pressure, :compute_volumetric_strain_rate)
    @eval begin
        @inline generate_equations(::AbstractCompositeModel, ::typeof($fn), ::Integer, ::Val, ::Val{false}, el_num = nothing; kwargs...) = ()
    end
end
@inline generate_equations(::Tuple{}; kwargs...) = ()

####

fn_pairs = (
    (compute_strain_rate, compute_stress),
    (compute_volumetric_strain_rate, compute_pressure),
    (compute_lambda, compute_plastic_strain_rate),
)

for pair in fn_pairs
    @eval begin
        @inline counterpart(::typeof($(pair[1]))) = $(pair[2])
        @inline counterpart(::typeof($(pair[2]))) = $(pair[1])
    end
end

@generated function get_own_functions(c::NTuple{N, AbstractCompositeModel}) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> get_own_functions(c[i])
    end
end
@inline get_own_functions(c::SeriesModel) = get_own_functions(c, series_state_functions, global_series_state_functions, local_series_state_functions)
@inline get_own_functions(c::ParallelModel) = get_own_functions(c, parallel_state_functions, global_parallel_state_functions, local_parallel_state_functions)

@inline function get_own_functions(c::AbstractCompositeModel, fn_state::F1, fn_global::F2, fn_local::F3) where {F1, F2, F3}
    fns_own_all = fn_state(c.leafs)
    fns_own_global = fn_global(fns_own_all) |> superflatten |> flatten_repeated_functions
    fns_own_local = fn_local(fns_own_all)
    return fns_own_global, fns_own_local
end

get_own_functions(::Tuple{}) = (), ()

@inline global_eltype_numbering(c::AbstractCompositeModel) = global_eltype_numbering(c, Ref(0), Ref(0), Ref(0))

@inline function global_eltype_numbering(c::AbstractCompositeModel, counter_v::Base.RefValue, counter_el::Base.RefValue, counter_pl::Base.RefValue)
    n1 = global_eltype_numbering(c.leafs, counter_v, counter_el, counter_pl)
    n2 = global_eltype_numbering(c.branches, counter_v, counter_el, counter_pl)
    return (n1, n2)
end

@generated function global_eltype_numbering(c::NTuple{N, AbstractRheology}, counter_v::Base.RefValue, counter_el::Base.RefValue, counter_pl::Base.RefValue) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> global_eltype_numbering(c[i], counter_v, counter_el, counter_pl)
    end
end

@generated function global_eltype_numbering(c::NTuple{N, AbstractCompositeModel}, counter_v::Base.RefValue, counter_el::Base.RefValue, counter_pl::Base.RefValue) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> global_eltype_numbering(c[i], counter_v, counter_el, counter_pl)
    end
end

@inline global_eltype_numbering(::Tuple{}, ::Base.RefValue, ::Base.RefValue, ::Base.RefValue) = ()

global_eltype_numbering(c::AbstractViscosity, counter_v::Base.RefValue, counter_el::Base.RefValue, counter_pl::Base.RefValue) = counter_v[] += 1
global_eltype_numbering(c::AbstractElasticity, counter_v::Base.RefValue, counter_el::Base.RefValue, counter_pl::Base.RefValue) = counter_el[] += 1
global_eltype_numbering(c::AbstractPlasticity, counter_v::Base.RefValue, counter_el::Base.RefValue, counter_pl::Base.RefValue) = counter_pl[] += 1

function add_global_equations(iparent, ilocal_childs, iself_ref, fns_own_global::F, leafs, ind_input, ::Val{B}, ::Val{1}, el_number) where {F, B}
    @inline
    iself_ref[] += 1
    return CompositeEquation(iparent, ilocal_childs, iself_ref[], fns_own_global, leafs, ind_input, Val(B), el_number)
end

@generated function add_local_equations(iparent, ilocal_childs, iself_ref, fns_own_local::F1, fns_own_global::F2, leafs, ::Val{N}, el_number) where {N, F1, F2}
    return quote
        Base.@ntuple $N i -> begin
            @inline
            iself_ref[] += 1
            add_local_equation(iparent, ilocal_childs, iself_ref[], fns_own_local[i], fns_own_global, leafs, 0, Val(false), el_number)
        end
    end
end

add_local_equations(::Any, ::Any, ::Any, ::F, ::F, ::Any, ::Any, ::Any) where {F} = ()

function add_local_equation(iparent, ilocal_childs, iself_ref, fns_own_local::F1, ::F2, leafs, num, ::Val{B}, el_number) where {F1 <: Function, F2 <: Function, B}
    return CompositeEquation(iparent, ilocal_childs, iself_ref, fns_own_local, leafs, num, Val(B), el_number)
end

## Exceptions
# lambda has no children or parent equations
# function add_local_equation(iparent, ilocal_childs, iself_ref, fns_own_local::typeof(compute_lambda), ::F, leafs, num, ::Val{B}, el_number) where {F<:Function, B}
#     CompositeEquation(0, (), iself_ref, fns_own_local, leafs, num, Val(true), el_number)
# end

add_local_equation(::Any, ::Any, ::Any, ::typeof(compute_lambda), ::typeof(compute_volumetric_strain_rate), ::Any, ::Any, ::Val{B}, ::Any) where {B} = ()

"""
    generate_args_template(eqs)
    generate_args_template(eqs, x, others)

Build `NamedTuple` argument templates for each equation. With `x` and `others`,
the returned tuple contains the differentiable value for the equation itself
plus all other differentiable values and nondifferentiable auxiliary fields.
"""
@generated function generate_args_template(eqs::NTuple{N, CompositeEquation}) where {N}
    return quote
        args = Base.@ntuple $N i -> differentiable_kwargs(eqs[i].fn)
    end
end

@generated function generate_args_template(eqs::NTuple{N, Any}, x::SVector{N}, others::NamedTuple) where {N}
    return quote
        args_template = generate_args_template(eqs)
        args = Base.@ntuple $N i -> begin
            @inline
            name = keys(args_template[i])
            merge(NamedTuple{name}(x[i]), others)
        end

        Base.@ntuple $N i -> begin
            diffs = Base.@ntuple $N j -> begin
                Base.structdiff(args[j], args[i])
            end
            merge(args[i], diffs...)
        end

    end
end

"""
    extract_local_kwargs(others, keys_hist, n)

Return `others` with history fields indexed to element `n`. Fields whose names
appear in `keys_hist` are treated as element-local tuples; tuple-valued fields
not listed in `keys_hist` use their first entry.

# Example
```julia
others = (; dt = 1e10, τ0 = (1.1, 3.0), d = (4, 2))
extract_local_kwargs(others, (:τ0,), 2) # τ0 = 3.0, d = 4
extract_local_kwargs(others, (:d,), 2)  # τ0 = 1.1, d = 2
```
"""
function extract_local_kwargs(others::NamedTuple, keys_hist::NTuple{M, Symbol}, n::Int64) where {M}
    vals_new = extract_local_kwargs(keys(others), values(others), keys_hist, n)
    return NamedTuple{keys(others)}(vals_new)
end

@generated function extract_local_kwargs(keys_args::NTuple{N, Symbol}, vals_args::NTuple{N, Union{_T, Tuple}}, keys_hist::NTuple{M, Symbol}, n::Int64) where {N, M, _T}
    return quote
        @inline
        Base.@ntuple $N i -> @inbounds _extract_local_kwargs(vals_args[i], keys_args[i], keys_hist, n)
    end
end

Base.@propagate_inbounds @inline _extract_local_kwargs(vals_args::Tuple, name, keys_hist, n) = ismember(name, keys_hist) ? vals_args[n] : vals_args[1]

@inline _extract_local_kwargs(vals_args, ::Any, ::Any, ::Any) = vals_args

@inline ismember(name::Symbol, keys_hist::NTuple{N, Symbol}) where {N} = name in keys_hist

"""
    evaluate_state_function(eq, args, others)
    evaluate_state_functions(eqs, args, others)

Evaluate one equation, or all equations, by calling each equation's state
function on its rheology elements and summing the element contributions.
Element-local history fields are extracted from `others` before dispatch.
"""
@generated function evaluate_state_functions(eqs::NTuple{N, CompositeEquation}, args, others) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> evaluate_state_function(eqs[i], args[i], others)
    end
end

"""
    evaluate_state_function(eq, args, others)

Evaluate one equation by calling its state function on each rheology element and
summing the element contributions. Element-local history fields are extracted
from `others` before dispatch.
"""
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

@inline evaluate_state_function(fn::F, rheology::Tuple{}, args, others) where {F} = 0.0e0

@generated function add_children(residual::NTuple{N, Any}, x::SVector{N}, eqs::NTuple{N, CompositeEquation}) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> residual[i] + add_child(x, eqs, eqs[i].child)
    end
end

@generated function add_child(x::SVector{M}, eqs::NTuple{N, CompositeEquation}, child::NTuple{NC}) where {M, N, NC}
    return quote
        @inline
        v = Base.@ntuple $NC i -> begin
            eq_ind = child[i]
            add_child(x, eqs[eq_ind], eq_ind)
        end
        sum(v)
    end
end

add_child(x, ::CompositeEquation, eq_ind) = x[eq_ind]
add_child(::SVector{N, T}, ::CompositeEquation{A, B, typeof(compute_lambda)}, eq_ind) where {N, A, B, T} = zero(T)
add_child(::SVector{N, T}, ::CompositeEquation{A, B, typeof(compute_lambda_parallel)}, eq_ind) where {N, A, B, T} = zero(T)
add_child(::SVector{N, T}, ::CompositeEquation{A, B, typeof(compute_plastic_strain_rate)}, eq_ind) where {N, A, B, T} = zero(T)
# add_child(::SVector{Any, T}, ::CompositeEquation{Any, Any, typeof(compute_plastic_strain_rate)}, eq_ind) where T = zero(T)

add_child(::SVector, ::Tuple{}) = 0.0e0
add_child(::SVector, ::NTuple{N, CompositeEquation}, ::Tuple{}) where {N} = 0.0e0

# if global, subtract the variables
@generated function subtract_parent(residual::NTuple{N, Any}, x, eqs::NTuple{N, CompositeEquation}, vars) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> begin
            residual[i] - subtract_parent(x, eqs[i], vars)
        end
    end
end

@inline subtract_parent(::SVector, eq::CompositeEquation{true}, vars) = vars[eq.ind_input]
@inline subtract_parent(x::SVector, eq::CompositeEquation{false}, ::NamedTuple) = x[eq.parent]
# exception for lambda
@inline subtract_parent(x::SVector, eq::CompositeEquation{false, T, typeof(compute_lambda)}, ::NamedTuple) where {T} = 0 # x[eq.self]
@inline subtract_parent(x::SVector, eq::CompositeEquation{false, T, typeof(compute_lambda_parallel)}, ::NamedTuple) where {T} = 0 # x[eq.self]

subtract_parent(residual::Number, x::SVector, eq::CompositeEquation, vars) = residual - subtract_parent(x, eq, vars)

"""
    compute_residual(c, x, vars, others)

Evaluate the residual vector for composite model `c` at solver vector `x`.
`vars` contains prescribed inputs such as strain-rate invariants; `others`
contains nondifferentiable auxiliary fields such as `dt`, `τ0`, `P0`, or other
history/state parameters.
"""
function compute_residual(c, x::SVector{N, T}, vars, others) where {N, T}

    eqs = generate_equations(c)
    @assert length(eqs) == length(x)
    args_all = generate_args_template(eqs, x, others)

    # evaluates the self-components of the residual
    residual = evaluate_state_functions(eqs, args_all, others)
    residual = add_children(residual, x, eqs)
    residual = subtract_parent(residual, x, eqs, vars)

    return SA[residual...]
end

# SeriesModel specialisation: subtracts the implicit elastic backstress correction
# from the global equation residual. For non-elastic composites the fallback in
# subtract_elastic_correction is a no-op, so there is no overhead.
function compute_residual(c::SeriesModel, x::SVector{N, T}, vars, others) where {N, T}

    eqs = generate_equations(c)
    @assert length(eqs) == length(x)
    args_all = generate_args_template(eqs, x, others)

    residual = evaluate_state_functions(eqs, args_all, others)
    residual = add_children(residual, x, eqs)
    residual = subtract_parent(residual, x, eqs, vars)
    residual = subtract_elastic_correction(c, eqs, residual, x, others)

    return SA[residual...]
end

# Subtract the implicit elastic correction from the first (global) residual entry.
# The correction is a scalar function of x (through the branch strain rates), so
# ForwardDiff differentiates through it automatically.
@inline function subtract_elastic_correction(c::SeriesModel, eqs, residual::NTuple{N}, x::SVector{N}, others) where {N}
    iselastic(c) == Val(false) && return residual
    cor = _implicit_elastic_correction(c, eqs, x, others)
    return _subtract_first(residual, cor)
end

# Return a new NTuple with the first entry decreased by `cor`.
@inline _subtract_first(r::NTuple{N}, cor) where {N} = (r[1] - cor, Base.tail(r)...)

function compute_residual(c, x::SVector{N, T}, vars, others, ::Int64, ::Int64) where {N, T}

    eqs = generate_equations(c)
    args_all = first(generate_args_template(eqs, x, others))

    # evaluates the self-components of the residual
    eq = first(eqs)
    residual = evaluate_state_function(eq, args_all, others)
    residual = add_children(residual, x, eq)
    residual = subtract_parent(residual, x, eq, vars)

    return residual
end
