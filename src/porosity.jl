
export compute_porosity

"""
    compute_porosity(c::SeriesModel, args)
    compute_porosity(elements::NTuple{N, AbstractRheology}, args)
    compute_porosity(r, args)

Return the updated porosity for a rheology model from element-wise porosity-rate
contributions.

`args` is a `NamedTuple` containing the local state needed by the participating
rheologies. At minimum, the tuple must contain:

- `Φ0`: previous porosity.
- `dt`: time step.

Additional fields, such as `p̄`, `pf`, previous pressures, or `λ`, are forwarded
to each element's `compute_porosity_rate` method. Elements that do not define a
porosity-rate contribution return zero through the fallback method.
"""
@inline compute_porosity(r::SeriesModel, args) = compute_porosity(r.leafs, args)

"""
    compute_porosity(r, args)

Fallback for rheologies without porosity evolution.
"""
@inline compute_porosity(r::T, args) where T = 0

"""
    compute_porosity(elements::NTuple{N, AbstractRheology}, args)

Sum the porosity-rate contributions of all rheologies in `elements` and advance
porosity explicitly as `Φ = Φ0 + dΦdt * dt`.
"""
@generated function compute_porosity(branches::NTuple{N, AbstractRheology}, args) where {N}
    quote
        # Keep the per-element calls unrolled so the generated code remains
        # type-stable for heterogeneous rheology tuples.
        dΦᵢdt = Base.@ntuple $N i -> compute_porosity_rate(branches[i]; args...)
        return args.Φ0 + sum(dΦᵢdt) * args.dt
    end
end

