using StaticArrays

const εxx_pure_shear = (1.0, -1.0, 0.0)

second_invariant_2D(ε) = sqrt((ε[1]^2 + ε[2]^2) / 2 + ε[3]^2)

function tensor_strain_rate_2D(εII; direction = εxx_pure_shear)
    directionII = second_invariant_2D(direction)
    return @. εII * direction / directionII
end

vars_2D(εII, θ = 0.0; direction = εxx_pure_shear) = (; ε = tensor_strain_rate_2D(εII; direction), θ)

zero_stress_tensor_2D() = (0.0, 0.0, 0.0)

function stress_tensor_from_invariant_2D(τII, ε_eff)
    εII = second_invariant_2D(ε_eff)
    η_eff = iszero(εII) ? zero(τII) : τII / (2 * εII)
    return Tuple(@. 2 * η_eff * ε_eff)
end

function elastic_stress_history_2D(c, τII, ε, τ0, others)
    ε_corr = effective_strain_rate_correction(c, ε, τ0, others)
    ε_eff  = ε .+ ε_corr
    return (stress_tensor_from_invariant_2D(τII, ε_eff),)
end

function elastic_stress_history_2D(c, x::SVector, ε, τ0, others)
    ε_corr = effective_strain_rate_correction(c, ε, τ0, others)
    ε_eff  = ε .+ ε_corr
    τII_elastic = compute_stress_elastic(c, x, others)
    return ntuple(i -> stress_tensor_from_invariant_2D(τII_elastic[i], ε_eff), length(τII_elastic))
end
