
@inline effective_strain_rate_correction(c::SeriesModel, vars, others) = effective_strain_rate_correction(c.leafs, c.branches, vars, others)

@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, vars, others) where {N}
    quote
        i = 0
        ε_elastic_cor = zero(vars.ε)
        Base.@nexprs $N j -> begin
            if isa(leafs[j], AbstractElasticity)
                i += 1
                # compute local effect on strainrate tensor due to old elastic stresses
                η = compute_viscosity_series(leafs[j], merge(vars, others))
                ε_elastic_cor += others.τ0[i] / (2 * η)
            end
        end
    
        return ε_elastic_cor
    end
end