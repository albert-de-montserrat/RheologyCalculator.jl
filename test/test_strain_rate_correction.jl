import RheologyCalculator:
    count_elastic, _n_elastic_in_parallel,
    _iselastic,
    _η_eff_maxwell, _η_eff_elastic, _η_KV,
    _weighted_backstress,
    _kv_branch_correction, _kv_corrections,
    effective_strain_rate_correction

# -----------------------------------------------------------------------
# Shared test parameters
# -----------------------------------------------------------------------
const η_v1 = 1e22
const η_v2 = 1e21
const η_v3 = 5e20
const G1   = 30e9
const G2   = 20e9
const dt   = 1e10
const εII  = 1e-14
const τ0_a = 2e6   # first backstress
const τ0_b = 1e6   # second backstress

# Element instances reused across tests
const visc1 = LinearViscosity(η_v1)
const visc2 = LinearViscosity(η_v2)
const visc3 = LinearViscosity(η_v3)
const el1   = IncompressibleElasticity(G1)
const el2   = IncompressibleElasticity(G2)

# args NamedTuple expected by internal viscosity dispatch
const args_base = (; ε = εII, dt = dt)

# -----------------------------------------------------------------------
# count_elastic
# -----------------------------------------------------------------------
@testset "count_elastic" begin
    # No elastic elements
    @test count_elastic((visc1,))           == 0
    @test count_elastic((visc1, visc2))     == 0
    @test count_elastic(())                 == 0

    # One elastic element
    @test count_elastic((el1,))             == 1
    @test count_elastic((visc1, el1))       == 1
    @test count_elastic((el1, visc1))       == 1

    # Two elastic elements
    @test count_elastic((el1, el2))         == 2
    @test count_elastic((visc1, el1, el2))  == 2
end

# -----------------------------------------------------------------------
# _n_elastic_in_parallel
# -----------------------------------------------------------------------
@testset "_n_elastic_in_parallel" begin
    # One direct elastic leaf, no sub-branches
    p1 = ParallelModel(visc2, el1)
    @test _n_elastic_in_parallel(typeof(p1)) == 1

    # No elastic elements at all
    p2 = ParallelModel(visc1, visc2)
    @test _n_elastic_in_parallel(typeof(p2)) == 0

    # Elastic only inside a Maxwell sub-branch
    p3 = ParallelModel(visc2, SeriesModel(visc3, el1))
    @test _n_elastic_in_parallel(typeof(p3)) == 1

    # Direct elastic leaf AND elastic in a sub-branch
    p4 = ParallelModel(el1, SeriesModel(visc3, el2))
    @test _n_elastic_in_parallel(typeof(p4)) == 2

    # Two direct elastic leafs
    p5 = ParallelModel(el1, el2)
    @test _n_elastic_in_parallel(typeof(p5)) == 2
end

# -----------------------------------------------------------------------
# _iselastic for composite model tuples
# -----------------------------------------------------------------------
@testset "_iselastic composite tuples" begin
    # Parallel branch that contains an elastic element
    p_elastic    = ParallelModel(visc2, el1)
    p_no_elastic = ParallelModel(visc1, visc2)

    @test _iselastic((p_elastic,))             == true
    @test _iselastic((p_no_elastic,))          == false
    @test _iselastic((p_elastic, p_no_elastic)) == true
    @test _iselastic((p_no_elastic, p_elastic)) == true
    @test _iselastic((p_no_elastic, p_no_elastic)) == false

    # Tiebreaker: empty tuple is never elastic
    @test _iselastic(()) == false
end

# -----------------------------------------------------------------------
# _η_eff_maxwell  —  harmonic mean of leaf viscosities
# -----------------------------------------------------------------------
@testset "_η_eff_maxwell" begin
    # Two-element Maxwell branch: visc + elastic.
    # η_eff_M = 1 / (1/η + 1/(G*dt))
    η_el = G1 * dt
    expected = 1 / (1/η_v2 + 1/η_el)
    @test _η_eff_maxwell((visc2, el1), args_base) ≈ expected  rtol=1e-12

    # Three-element: visc1 + visc2 + elastic  (all in series inside a branch)
    expected3 = 1 / (1/η_v1 + 1/η_v2 + 1/(G1*dt))
    @test _η_eff_maxwell((visc1, visc2, el1), args_base) ≈ expected3  rtol=1e-12

    # Purely viscous tuple (no elastic leaf): η_eff_M = harmonic mean of two viscosities
    expected_v = 1 / (1/η_v1 + 1/η_v2)
    @test _η_eff_maxwell((visc1, visc2), args_base) ≈ expected_v  rtol=1e-12
end

# -----------------------------------------------------------------------
# _η_eff_elastic  —  G * dt for the elastic leaf
# -----------------------------------------------------------------------
@testset "_η_eff_elastic" begin
    # Elastic leaf is present: should return G*dt
    @test _η_eff_elastic((visc2, el1), args_base) ≈ G1 * dt  rtol=1e-12
    @test _η_eff_elastic((el1, visc2), args_base) ≈ G1 * dt  rtol=1e-12

    # No elastic leaf: should return 0.0
    @test _η_eff_elastic((visc1, visc2), args_base) == 0.0
end

# -----------------------------------------------------------------------
# _η_KV  —  arithmetic sum of effective viscosities
# -----------------------------------------------------------------------
@testset "_η_KV" begin
    η_el1 = G1 * dt
    η_el2 = G2 * dt

    # Leafs only, no sub-branches:  η_KV = η_v2 + G1*dt
    @test _η_KV((visc2, el1), (), args_base) ≈ η_v2 + η_el1  rtol=1e-12

    # Leafs only, two viscous:  η_KV = η_v1 + η_v2
    @test _η_KV((visc1, visc2), (), args_base) ≈ η_v1 + η_v2  rtol=1e-12

    # One viscous leaf + one Maxwell sub-branch:
    # η_KV = η_v2 + η_eff_M(visc3, el1)
    η_eff_M = 1 / (1/η_v3 + 1/η_el1)
    sub1 = SeriesModel(visc3, el1)
    @test _η_KV((visc2,), (sub1,), args_base) ≈ η_v2 + η_eff_M  rtol=1e-12

    # Two Maxwell sub-branches (no leaf viscosity), both use visc3:
    # sub1 = SeriesModel(visc3, el1), sub2 = SeriesModel(visc3, el2)
    # η_KV = η_eff_M(visc3, el1) + η_eff_M(visc3, el2)
    η_eff_M1 = 1 / (1/η_v3 + 1/η_el1)
    η_eff_M2 = 1 / (1/η_v3 + 1/η_el2)
    sub2 = SeriesModel(visc3, el2)
    @test _η_KV((), (sub1, sub2), args_base) ≈ η_eff_M1 + η_eff_M2  rtol=1e-12
end

# -----------------------------------------------------------------------
# _weighted_backstress
# -----------------------------------------------------------------------
@testset "_weighted_backstress" begin
    η_el1 = G1 * dt
    η_el2 = G2 * dt

    # Direct elastic leaf: η_star = 1, so ws = τ0[el_idx_start + 1]
    τ0 = (τ0_a,)
    @test only(_weighted_backstress((visc2, el1), (), εII, τ0, args_base, 0)) ≈ τ0_a  rtol=1e-12

    # No elastic elements: ws = 0
    τ0_empty = (τ0_a,)
    @test only(_weighted_backstress((visc1, visc2), (), εII, τ0_empty, args_base, 0)) ≈ 0.0

    # Maxwell sub-branch: η_star = η_eff_M / η_el
    sub1 = SeriesModel(visc3, el1)
    η_eff_M = 1 / (1/η_v3 + 1/η_el1)
    η_star  = η_eff_M / η_el1
    τ0_sub = (τ0_a,)
    @test only(_weighted_backstress((visc2,), (sub1,), εII, τ0_sub, args_base, 0)) ≈ η_star * τ0_a  rtol=1e-12

    # Direct elastic leaf + Maxwell sub-branch in same parallel block.
    # el_idx_start = 0: el1 → τ0[1]=τ0_a (η_star=1), sub1.el → τ0[2]=τ0_b (η_star=η_eff_M/η_el)
    sub1b = SeriesModel(visc3, el2)
    η_eff_M2 = 1 / (1/η_v3 + 1/η_el2)
    η_star2  = η_eff_M2 / η_el2
    τ0_two = (τ0_a, τ0_b)
    ws_expected = τ0_a + η_star2 * τ0_b
    @test only(_weighted_backstress((el1,), (sub1b,), εII, τ0_two, args_base, 0)) ≈ ws_expected  rtol=1e-12

    # Non-zero el_idx_start (a leaf elastic already consumed offset=1):
    # τ0[1] belongs to a previous leaf; this branch uses τ0[2]
    τ0_offset = (0.0, τ0_a)
    @test only(_weighted_backstress((visc2, el1), (), εII, τ0_offset, args_base, 1)) ≈ τ0_a  rtol=1e-12
end

# -----------------------------------------------------------------------
# effective_strain_rate_correction — end-to-end
# -----------------------------------------------------------------------
@testset "effective_strain_rate_correction — Maxwell leaf" begin
    # SeriesModel(visc1, el1): classic single-spring correction = τ0 / (2*G*dt)
    c = SeriesModel(visc1, el1)
    τ0     = (τ0_a,)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0,))
    expected = τ0_a / (2 * G1 * dt)
    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) ≈ expected  rtol=1e-12

    # No elastic elements: zero correction
    c2     = SeriesModel(visc1, visc2)
    others2 = (; dt = dt, τ0 = (0.0,), P0 = (0.0,))
    @test effective_strain_rate_correction(c2, ε, (0.0,), others2) == 0
end

@testset "effective_strain_rate_correction — KV branch (direct elastic leaf)" begin
    # SeriesModel(visc1, ParallelModel(visc2, el1))
    # η_KV = η_v2 + G1*dt,  correction = τ0 / (2*η_KV)   [η_star=1 for direct leaf]
    c = SeriesModel(visc1, ParallelModel(visc2, el1))
    τ0     = (τ0_a,)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0,))
    η_KV   = η_v2 + G1 * dt
    expected = τ0_a / (2 * η_KV)
    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) ≈ expected  rtol=1e-12
end

@testset "effective_strain_rate_correction — generalized Maxwell (Maxwell sub-branch)" begin
    # SeriesModel(visc1, ParallelModel(visc2, SeriesModel(visc3, el1)))
    # η_eff_M = 1/(1/η_v3 + 1/(G1*dt))
    # η_KV    = η_v2 + η_eff_M
    # η_star  = η_eff_M / (G1*dt)
    # correction = η_star * τ0 / (2*η_KV)
    c = SeriesModel(visc1, ParallelModel(visc2, SeriesModel(visc3, el1)))
    τ0     = (τ0_a,)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0,))
    η_el    = G1 * dt
    η_eff_M = 1 / (1/η_v3 + 1/η_el)
    η_KV    = η_v2 + η_eff_M
    η_star  = η_eff_M / η_el
    expected = η_star * τ0_a / (2 * η_KV)
    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) ≈ expected  rtol=1e-12
end

@testset "effective_strain_rate_correction — Maxwell leaf + generalized Maxwell branch" begin
    # SeriesModel(el1, ParallelModel(visc2, SeriesModel(visc3, el2)))
    # Leaf correction for el1:      τ0[1] / (2*G1*dt)
    # Branch correction for el2:    η_star * τ0[2] / (2*η_KV)
    c = SeriesModel(el1, ParallelModel(visc2, SeriesModel(visc3, el2)))
    τ0     = (τ0_a, τ0_b)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0, 0.0))
    η_el2   = G2 * dt
    η_eff_M = 1 / (1/η_v3 + 1/η_el2)
    η_KV    = η_v2 + η_eff_M
    η_star  = η_eff_M / η_el2
    leaf_cor   = τ0_a / (2 * G1 * dt)
    branch_cor = η_star * τ0_b / (2 * η_KV)
    expected   = leaf_cor + branch_cor
    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) ≈ expected  rtol=1e-12
end

@testset "effective_strain_rate_correction — two parallel branches" begin
    # SeriesModel(visc1,
    #   ParallelModel(visc2, el1),          ← branch 1: direct leaf, τ0[1]
    #   ParallelModel(visc3, SeriesModel(visc2, el2))  ← branch 2: Maxwell, τ0[2]
    # )
    c = SeriesModel(visc1,
        ParallelModel(visc2, el1),
        ParallelModel(visc3, SeriesModel(visc2, el2)))
    τ0     = (τ0_a, τ0_b)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0, 0.0))

    η_KV1    = η_v2 + G1 * dt
    cor1     = τ0_a / (2 * η_KV1)

    η_el2    = G2 * dt
    η_eff_M2 = 1 / (1/η_v2 + 1/η_el2)
    η_KV2    = η_v3 + η_eff_M2
    η_star2  = η_eff_M2 / η_el2
    cor2     = η_star2 * τ0_b / (2 * η_KV2)

    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) ≈ cor1 + cor2  rtol=1e-12
end

@testset "effective_strain_rate_correction — two elastic sources in one branch" begin
    # SeriesModel(visc1, ParallelModel(el1, SeriesModel(visc3, el2)))
    # Branch has one direct elastic leaf (el1) and one Maxwell sub-branch (visc3+el2).
    # η_KV   = G1*dt + η_eff_M(visc3, el2)
    # ws     = τ0[1]*1 + τ0[2]*η_star   (η_star = η_eff_M / (G2*dt))
    # correction = ws / (2*η_KV)
    c = SeriesModel(visc1, ParallelModel(el1, SeriesModel(visc3, el2)))
    τ0     = (τ0_a, τ0_b)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0, 0.0))

    η_el2   = G2 * dt
    η_eff_M = 1 / (1/η_v3 + 1/η_el2)
    η_KV    = G1 * dt + η_eff_M
    η_star  = η_eff_M / η_el2
    expected = (τ0_a + η_star * τ0_b) / (2 * η_KV)
    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) ≈ expected  rtol=1e-12
end

@testset "effective_strain_rate_correction — zero backstress gives zero correction" begin
    c      = SeriesModel(visc1, ParallelModel(visc2, SeriesModel(visc3, el1)))
    τ0     = (0.0,)
    ε      = (εII,)
    others = (; dt = dt, τ0 = τ0, P0 = (0.0,))
    cor = effective_strain_rate_correction(c, ε, τ0, others)
    @test only(cor) == 0.0
end
