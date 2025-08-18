@inline function augment_args(args, Δx)
    k = keys(args)
    vals = MVector(values(args))
    for i in eachindex(Δx)
        vals[i] += Δx[i]
    end
    return (; zip(k, vals)...)
end

@inline function update_args2(args, x::SVector{N, T}) where {N, T}
    k = keys(args)
    N0 = length(args)
    vals = @MVector zeros(T, N0)
    for i in 1:length(args)
        vals[i] = x[i]
    end
    return (; zip(k, vals)...)
end

# Define history kwargs variables (which we can define as tuples )
history_kwargs(::AbstractElasticity) = (:τ0, :P0)
history_kwargs(::AbstractViscosity) = (:d,)
history_kwargs(::AbstractPlasticity) = ()

# dummy NamedTuple allocators
@inline residual_kwargs(::Type{T}, ::Function) where {T} = (; tmp = zero(T))
@inline residual_kwargs(::Type{T}, ::typeof(compute_strain_rate)) where {T} = (; ε = zero(T))
@inline residual_kwargs(::Type{T}, ::typeof(compute_volumetric_strain_rate)) where {T} = (; θ = zero(T))
@inline residual_kwargs(::Type{T}, ::typeof(compute_stress)) where {T} = (; τ = zero(T))
@inline residual_kwargs(::Type{T}, ::typeof(compute_pressure)) where {T} = (; P = zero(T))

@inline residual_kwargs(funs::F) where {F <: Function} = residual_kwargs(Float64, funs)
@inline residual_kwargs(funs::NTuple{N, Any}) where {N} = residual_kwargs.(Float64, funs)

# dummy NamedTuple allocators
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_strain_rate)) where {T} = (; τ = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_volumetric_strain_rate)) where {T} = (; P = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_lambda)) where {T} = (; λ = zero(T)) # τ = zero(T), P = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_stress)) where {T} = (; ε = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_pressure)) where {T} = (; θ = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_plastic_strain_rate)) where {T} = (; τ_pl = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_plastic_stress)) where {T} = (; τ_pl = zero(T))
@inline differentiable_kwargs(::Type{T}, ::typeof(compute_volumetric_plastic_strain_rate)) where {T} = (; τ_pl = zero(T), P_pl = zero(T))
@inline differentiable_kwargs(fun::F) where {F <: Function} = differentiable_kwargs(Float64, fun)

@inline differentiable_kwargs(::Tuple{}) = (;)
@inline differentiable_kwargs(funs::NTuple{N, Any}) where {N} = differentiable_kwargs(Float64, funs)

@generated function differentiable_kwargs(::Type{T}, funs::NTuple{N, Any}) where {N, T}
    return quote
        @inline
        Base.@nexprs $N i -> nt_i = differentiable_kwargs($T, funs[i])
        Base.@ncall $N merge nt
    end
end
differentiable_kwargs(::Type{T}, funs::NTuple{0, Any}) where {T} = (;)

@inline all_differentiable_kwargs(funs::NTuple{N, Any}) where {N} = all_differentiable_kwargs(Float64, funs)

@generated function all_differentiable_kwargs(::Type{T}, funs::NTuple{N, Any}) where {N, T}
    return quote
        @inline
        Base.@ntuple $N i -> differentiable_kwargs($T, funs[i])
    end
end

function split_args(args, statefuns::NTuple{N, Any}) where {N}
    # split args into differentiable and not differentiable
    dummy = differentiable_kwargs(statefuns)
    args_nondiff = Base.structdiff(args, dummy)
    args_diff = Base.structdiff(dummy, args_nondiff)
    args_diff0 = Base.structdiff(args, args_nondiff)

    args_diff = merge(args_diff, args_diff0)
    return args_diff, args_nondiff
end
