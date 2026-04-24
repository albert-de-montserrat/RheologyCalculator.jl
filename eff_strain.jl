# 2 springs, one dashpot, one parallel spring/dashpot
exx    = -1e-14
exy    = 0.5e-15

txx10   =-1e6
txy10   = 0.5e6

txx20   =-2e6
txy20   = 0.8e6

G1     = 1e10
G2     = 2e10
G3     = 3e10
dt     = 1e10
eta_v1  = 1e20
eta_v2  = 4e20

eta_e1  = G1*dt
eta_e2  = G2*dt
eta_e3  = G3*dt

eta_K  = eta_v2 + eta_e3
eta_ve  = 1/(1/eta_e1 + 1/eta_e2 + 1/eta_v1 + 1/eta_K)

txx    = 2*eta_ve*(exx + txx10/2/eta_e1 + txx10/2/eta_e2 + txx20/2/eta_K)
txy    = 2*eta_ve*(exy + txy10/2/eta_e1 + txy10/2/eta_e2 + txy20/2/eta_K)
tII    = sqrt(txx^2 + txy^2)

eII    = sqrt(exx^2 + exy^2)

exx_v = txx/2/eta_v1
exy_v = txy/2/eta_v1
eII_v = sqrt(exx_v^2 + exy_v^2)

exx_e1 = (txx-txx10)/2/eta_e1
exy_e1 = (txy-txy10)/2/eta_e1
eII_e1 = sqrt(exx_e1^2 + exy_e1^2)

exx_e2 = (txx-txx10)/2/eta_e2
exy_e2 = (txy-txy10)/2/eta_e2
eII_e2 = sqrt(exx_e2^2 + exy_e2^2)

exx_K  = (txx-txx20)/2/eta_K
exy_K  = (txy-txy20)/2/eta_K
eII_K  = sqrt(exx_K^2 + exy_K^2)

display((exx - exx_v - exx_e1 - exx_e2 - exx_K)/eII)
display((exy - exy_v - exy_e1 - exy_e2 - exy_K)/eII)
display((eII - eII_v - eII_e1 - eII_e2 - eII_K)/eII)

exx_eff = exx + txx10/2/eta_e1 + txx10/2/eta_e2 + txx20/2/eta_K
exy_eff = exy + txy10/2/eta_e1 + txy10/2/eta_e2 + txy20/2/eta_K
eII_eff = sqrt(exx_eff^2 + exy_eff^2)
display((eII_eff - tII/2/eta_v1 - tII/2/eta_e1 - tII/2/eta_e2 - tII/2/eta_K)/eII)


@generated function effective_strain_rate_correction(c::SeriesModel, ε::NTuple{N}, τ0::NTuple{N}, others) where N
    quote
       Base.@ntuple $N i -> effective_strain_rate_correction(c, ε[i], τ0[i], others)
    end
end

@inline effective_strain_rate_correction(c::SeriesModel, ε, τ0, others) = effective_strain_rate_correction(c.leafs, c.branches, ε, τ0, others)

@generated function effective_strain_rate_correction(leafs::NTuple{N, Any}, ::Tuple{}, ε, τ0, others) where {N}
    quote
        i = 0
        ε_elastic_cor = zero(ε)
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