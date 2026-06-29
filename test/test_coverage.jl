import RheologyCalculator: compute_stress_elastic

# ── compute_stress_elastic ────────────────────────────────────────────────────
# compute_pressure_elastic is already exercised in test_ModCamClay / test_Hyperbolic /
# test_VEPCap. compute_stress_elastic is imported in runtests.jl but never called.

@testset "compute_stress_elastic" begin
    c      = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
    vars   = (; ε = 1e-15)
    args   = (; τ = 0.0)
    others = (; dt = 1e10, τ0 = (0.0,), P0 = (0.0,))
    x      = initial_guess_x(c, vars, args, others)
    x      = solve(c, x, vars, others)
    τ_el   = compute_stress_elastic(c, x, others)
    @test length(τ_el) == 1
    @test isfinite(only(τ_el))
    @test only(τ_el) ≈ x[1]   # strain-rate-form elastic: x[self] is the stress
end

# ── solve verbose path (it > 1) ───────────────────────────────────────────────
# The println("Iterations: ...") branch is gated by `verbose && it > 1`.
# PowerLaw is nonlinear: starting from x = SA[1e5] (≈ 1000x the equilibrium
# stress of ~100 Pa), the initial residual is ~1e-5 >> atol=1e-12, forcing
# multiple Newton iterations so the branch is reached.

@testset "solve verbose" begin
    c      = SeriesModel(PowerLawViscosity(5e19, 3))
    vars   = (; ε = 1e-14)
    others = (; τ0 = (0.0,))  # no elastic elements; τ0 must be a non-empty NTuple
    τ_eq   = (2 * 5e19 * 1e-14)^(1/3)  # ≈ 100 Pa
    xnorm  = normalisation_x(c, τ_eq, 1e-14)  # proper scale → relative tolerance kicks in
    x      = SA[1e5]                           # coarse start forces many Newton steps (it > 1)
    x      = solve(c, x, vars, others; verbose = true, xnorm0 = xnorm)
    @test x[1] ≈ τ_eq rtol = 1e-3
end
