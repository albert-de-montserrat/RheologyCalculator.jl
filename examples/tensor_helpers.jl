const εxx_pure_shear = (1.0, -1.0, 0.0)

second_invariant_2D(ε) = sqrt((ε[1]^2 + ε[2]^2) / 2 + ε[3]^2)

function tensor_strain_rate_2D(εII; direction = εxx_pure_shear)
    directionII = second_invariant_2D(direction)
    return @. εII * direction / directionII
end

zero_stress_tensor_2D() = (0.0, 0.0, 0.0)

function elastic_stress_history_2D(c, τII, ε, τ0, others)
    ε_corr = effective_strain_rate_correction(c, ε, τ0, others)
    ε_eff  = ε .+ ε_corr
    εII    = second_invariant_2D(ε_eff)
    η_eff  = iszero(εII) ? zero(τII) : τII / (2 * εII)
    return (Tuple(@. 2 * η_eff * ε_eff),)
end
