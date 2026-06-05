
export compute_strain_rate, compute_stress
export compute_viscosity_strain, compute_viscosity_stress

@inline function compute_viscosity_stress(r::LinearViscosityStress, τ; kwargs...) 
    ε = compute_strain_rate()
end

@generated function compute_viscosity_Kelvin(branches::NTuple{N, AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            # compute local effect on strainrate tensor due to old elastic stresses
            η = compute_viscosity_Kelvin(branches[j], merge((; ε, τ), others))
            η_eff += η
        end
        return η_eff
    end
end

@generated function compute_viscosity_Maxwell(branches::NTuple{N, AbstractRheology}, ε, τ, others) where {N}
    quote
        η_eff = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            # compute local effect on strainrate tensor due to old elastic stresses
            η = compute_viscosity_Maxwell(branches[j], merge((; ε, τ), others))
            η_eff += 1 / η
        end
        return η_eff
    end
end
