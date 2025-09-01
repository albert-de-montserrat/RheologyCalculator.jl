
# @inline effective_strain_rate_correction(c::SeriesModel, vars, others) = effective_strain_rate_correction(c.leafs, c.branches, vars, others)

# @generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, vars, others) where {N}
#     quote
#         i = 0
#         ε_elastic_cor = zero(vars.ε)
#         Base.@nexprs $N j -> begin
#             if isa(leafs[j], AbstractElasticity)
#                 i += 1
#                 # compute local effect on strainrate tensor due to old elastic stresses
#                 η = compute_viscosity_series(leafs[j], merge(vars, others))
#                 ε_elastic_cor += others.τ0[i] / (2 * η)
#             end
#         end
    
#         return ε_elastic_cor
#     end
# end

@generated function effective_strain_rate_correction(c::SeriesModel, ε::NTuple{N}, τ0::NTuple{N}, others) where N
    quote
       Base.@ntuple $N i -> effective_strain_rate_correction(c, ε[i], τ0[i], others)
    end
end

@inline effective_strain_rate_correction(c::SeriesModel, ε::Number, τ0, others) = effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)

@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0, others) where {N}
    quote
        i = 0
        ε_elastic_cor = zero(eltype(ε))
        Base.@nexprs $N j -> begin
            if isa(leafs[j], AbstractElasticity)
                i += 1
                # compute local effect on strainrate tensor due to old elastic stresses
                η = compute_viscosity_series(leafs[j], merge((; ε), others))
                ε_elastic_cor += τ0[i] / (2 * η)
            end
        end
    
        return ε_elastic_cor
    end
end

@inline second_invariant(xx, yy, xy) = √(0.5 * (xx^2 + yy^2) + xy^2)
@inline second_invariant(xx, yy, zz, yz, xz, xy) = √(0.5 * (xx^2 + yy^2 + zz^2) + xy^2 + yz^2 + xz^2)