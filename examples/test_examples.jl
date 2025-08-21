using Test, Statistics
# Runs all tests in the /examples directory

include("Burgers.jl")
@test  mean(τ1) ≈ 148614.80030508823
@test  mean(τ2) ≈ 77593.22577542417
@test  mean(P1) ≈ -2.057142856796177

include("Maxwell_VE.jl")
@test mean(abs.(τ-τ_an)) ≈ 9993.326538891743

include("Maxwell_VEP.jl")
@test mean(abs.(τ-τ_an)) ≈ 70.91382257194665

include("Maxwell_VEVP.jl")
@test mean(τ) ≈ 9.413876094263624e6

include("Maxwell_VEPCap.jl")
@test mean(τ1) ≈ 0.0
@test mean(P1) ≈ -348008.4518961142

@test mean(τ2) ≈ 543634.3867835554
@test mean(P2) ≈  33645.5280111008

@test mean(τ3) ≈ 397093.10397308524
@test mean(P3) ≈ -243552.39195758212

#include("transient_creep.jl")
##@test mean(τ1) ≈ 136871.94383758004
#@test mean(τ2) ≈ 103128.05616241995




