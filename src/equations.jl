# include("recursion.jl")

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
        #return new{B, T, F, R, RT}(parent, child, self, fn, rheology, ind_input, el_number)
        return new{B, T, F, R, RT}(parent, child, self, fn, rheology, ind_input, el_number)
        
    end
end

@inline generate_equations(c::AbstractCompositeModel) = generate_equations(c, global_series_functions(c))

@generated function generate_equations(c::AbstractCompositeModel, fns::NTuple{N, Any}) where {N}
    return quote
        @inline
        iparent = 0
        iself = 0
        isGlobal = Val(true)
        el_num = global_eltype_numbering(c) # global element numbering (to be followed )
        eqs = Base.@ntuple $N i -> begin
            ind_input = i
            eqs = generate_equations(c, fns[i], ind_input, isGlobal, isvolumetric(c), el_num; iparent = iparent, iself = iself)
            iself = eqs[end].self
            iparent = 0
            eqs
        end
        superflatten(eqs)
    end
end

function generate_equations(c::AbstractCompositeModel, fns_own_global::F, ind_input, ::Val{B}, ::Val, el_num = nothing; iparent::Int64 = 0, iself::Int64 = 0) where {F, B}
    iself_ref = Ref{Int64}(iself)
    (; branches, leafs) = c
    local_el = el_num[1]

    _, fns_own_local = get_own_functions(c)
    # fns_branches_global,_ = get_own_functions(branches)

    nown = 1 # length(fns_own_global)
    nlocal = length(fns_own_local)
    nbranches = length(branches)

    # iglobal          = ntuple(i -> iparent + i - 1, Val(nown))
    ilocal_childs = ntuple(i -> iself + nown + i, Val(nlocal))
    offsets_parallel = (0, ntuple(i -> i, Val(nbranches))...)
    # offsets_parallel = (0, length.(fns_branches_global)...)
    iparallel_childs = ntuple(i -> iself + nlocal + offsets_parallel[i] + i + nown, Val(nbranches))
    # ichildren = (ilocal_childs..., iparallel_childs...)

    # add globals
    # iself_ref[] += 1
    # global_eqs   = CompositeEquation(iparent, iparallel_childs, iself_ref[], fns_own_global, leafs, Val(false))
    isGlobal = Val(B)
    global_eqs0 = add_global_equations(iparent, ilocal_childs, iparallel_childs, iself_ref, fns_own_global, leafs, branches, ind_input, isGlobal, Val(1), local_el)
    local_eqs = add_local_equations(global_eqs0.self, (), iself_ref, fns_own_local, fns_own_global, leafs, Val(nlocal), local_el)
    # need to correct the children of the global equations in absence of local equations (i.e. remove them)
    global_eqs = correct_children(global_eqs0, local_eqs)

    fn = counterpart(fns_own_global)
    parallel_eqs = generate_equations_unroller(branches, fn, el_num, global_eqs, iself_ref)

    return (global_eqs, local_eqs..., parallel_eqs...) |> superflatten
end

@generated function generate_equations_unroller(branches::NTuple{N, Any}, fn::F, el_num, global_eqs, iself_ref) where {N, F}
    return quote
        @inline
        Base.@ntuple $N i -> begin
            generate_equations(branches[i], fn, 0, Val(false), isvolumetric(branches[i]), el_num[2][i]; iparent = global_eqs.self, iself = iself_ref[])
        end
    end
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

@inline get_own_functions(c::NTuple{N, AbstractCompositeModel}) where {N} = ntuple(i -> get_own_functions(c[i]), Val(N))
@inline get_own_functions(c::SeriesModel) = get_own_functions(c, series_state_functions, global_series_state_functions, local_series_state_functions)
@inline get_own_functions(c::ParallelModel) = get_own_functions(c, parallel_state_functions, global_parallel_state_functions, local_parallel_state_functions)

@inline function get_own_functions(c::AbstractCompositeModel, fn_state::F1, fn_global::F2, fn_local::F3) where {F1, F2, F3}
    fns_own_all = fn_state(c.leafs)
    fns_own_global = fn_global(fns_own_all) |> superflatten |> flatten_repeated_functions
    fns_own_local = fn_local(fns_own_all)
    return fns_own_global, fns_own_local
end

get_own_functions(::Tuple{}) = (), ()
# get_own_functions(::Tuple{}) = compute_strain_rate, ()

# # Number the rheological elements sequantially
# function global_el_numbering(c::NTuple{N, AbstractCompositeModel}, v=0) where N
#     # This allocates
#     n = ();
#     for i=1:N
#         nmax = maximum(superflatten(n), init=v)
#         loc_num = global_el_numbering(c[i], nmax)
#         n = (n..., loc_num)
#     end

#    return n
# end

# @generated function global_el_numbering(c::NTuple{N, AbstractRheology}, v=0) where N
#     #N = ntuple(i -> i + v, Val(N))
#     quote
#         Base.@ntuple $N i-> begin
#             @inline
#                 i + v
#         end
#     end
# end


# function global_el_numbering(c::AbstractCompositeModel,v=0)
#     num_leafs    = global_el_numbering(c.leafs,v)
#     num_branches = global_el_numbering(c.branches,maximum(num_leafs, init=v))
#     return num_leafs, num_branches
# end
# global_el_numbering(::Tuple{},v=0) = ()

@inline global_el_numbering(c::AbstractCompositeModel) = global_el_numbering(c, Ref(0))

@inline function global_el_numbering(c::AbstractCompositeModel, counter::Base.RefValue)
    n1 = global_el_numbering(c.leafs, counter)
    n2 = global_el_numbering(c.branches, counter)
    return (n1, n2)
end

@generated function global_el_numbering(::NTuple{N, AbstractRheology}, counter::Base.RefValue) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> counter[] += 1
    end
end

@generated function global_el_numbering(c::NTuple{N, AbstractCompositeModel}, counter::Base.RefValue) where {N}
    return quote
        @inline
        Base.@ntuple $N i -> global_el_numbering(c[i], counter)
    end
end

@inline global_el_numbering(::Tuple{}, ::Base.RefValue) = ()

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

#get_local_functions(c::NTuple{N, AbstractCompositeModel}) where N = ntuple(i -> get_own_functions(c[i]), Val(N))
function get_local_functions(c::SeriesModel)
    fns_own_all = series_state_functions(c.leafs)
    return local_series_state_functions(fns_own_all)
end

function get_local_functions(c::ParallelModel)
    fns_own_all = parallel_state_functions(c.leafs)
    return local_parallel_state_functions(fns_own_all)
end

@inline has_children(::F, branch) where {F} = Val(true)
@inline has_children(::typeof(compute_pressure), branch) = isvolumetric(branch)
@inline has_children(::typeof(compute_volumetric_strain_rate), branch) = isvolumetric(branch)

@inline correct_children(fn::F, branch::AbstractCompositeModel, children) where {F} = correct_children(children, has_children(fn, branch))
@generated function correct_children(fn::F, branch::NTuple{N, AbstractCompositeModel}, children) where {F, N}
    return quote
        new_children = Base.@ntuple $N i -> correct_children(children[i], has_children(fn, branch[i]))
        return superflatten(new_children)
    end
end
@inline correct_children(children, ::Val{true}) = children
@inline correct_children(::Any, ::Val{false}) = ()

@generated function add_global_equations(iparent, ilocal_childs, iparallel_childs, iself_ref, fns_own_global::NTuple{F, Any}, leafs, branches, ::Val{N}, el_number) where {F, N}
    return quote
        Base.@ntuple $N i -> begin
            @inline
            iself_ref[] += 1
            corrected_children = correct_children(fns_own_global[i], branches, iparallel_childs)
            children = (ilocal_childs..., corrected_children...) .+ (i - 1)
            CompositeEquation(iparent, children, iself_ref[], fns_own_global[i], leafs, Val(true), el_number)
        end
    end
end

@generated function add_global_equations(iparent, ilocal_childs, iparallel_childs, iself_ref, fns_own_global::F, leafs, branches, ::Val{N}, el_number) where {F, N}
    return quote
        Base.@ntuple $N i -> begin
            @inline
            iself_ref[] += 1
            corrected_children = correct_children(fns_own_global[i], branches, iparallel_childs)
            children = (ilocal_childs..., corrected_children...) .+ (i - 1)
            CompositeEquation(iparent, children, iself_ref[], fns_own_global, leafs, Val(true), el_number)
        end
    end
end

function add_global_equations(iparent, ilocal_childs, iparallel_childs, iself_ref, fns_own_global::F, leafs, branches, ind_input, ::Val{B}, ::Val{1}, el_number) where {F, B}
    @inline
    iself_ref[] += 1
    corrected_children = correct_children(fns_own_global, branches, iparallel_childs)
    children = (ilocal_childs..., corrected_children...)
    #@show children, corrected_children
    return CompositeEquation(iparent, children, iself_ref[], fns_own_global, leafs, ind_input, Val(B), el_number)
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

@generated function add_parallel_equations(global_eqs::NTuple{N1, Any}, branches::NTuple{N2, AbstractCompositeModel}, iself_ref, fns_own_global::NTuple{N1, Any}) where {N1, N2}
    return quote
        Base.@ntuple $N1 j -> begin
            @inline
            iparent_new = global_eqs[j].self
            fn = counterpart(fns_own_global[j])
            Base.@ntuple $N2 i -> begin
                @inline
                generate_equations(branches[i], fn, isvolumetric(branches[i]); iparent = iparent_new, iself = iself_ref[])
            end
        end
    end
end

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
        args_merged = merge(args...)
        Base.@ntuple $N i -> args_merged
    end
end

# extracts local kwargs when the args is a tuple, and it is listed
# julia> others      = (; dt = 1e10, τ0 = (1.1, 3.0), d = (4, 2))
# julia> extract_local_kwargs(others, (:τ0,), 2)
# (dt = 1.0e10, τ0 = 3.0, d = 4)
# julia> extract_local_kwargs(others, (:d,), 2)
# (dt = 1.0e10, τ0 = 1.1, d = 2)
# julia> extract_local_kwargs(others, (:τ0,:d), 2)
# (dt = 1.0e10, τ0 = 3.0, d = 2)
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
