using Test, Statistics
# Runs all tests in the /examples directory

include("Burgers.jl")
@test  mean(τ1) ≈ 148614.80030508823
@test  mean(τ2) ≈ 77593.22577542417
@test  mean(P1) ≈ -2.057142856796177


include("Maxwell_VE.jl")
@test mean(abs.(τ-τ_an)) ≈ 10022.097942015174

include("Maxwell_VEP.jl")
@test mean(abs.(τ-τ_an)) ≈ 70.92833227623228

include("Maxwell_VEVP.jl")
@test mean(τ) ≈ 1.8481717756725273e6

include("Mode1_Mode2_VEP.jl")
@test mean(τ1) ≈ 0.0
@test mean(P1) ≈ -348008.4518960728

@test mean(τ2) ≈ 543634.3867782074
@test mean(P2) ≈ 33645.52800041107

@test mean(τ3) ≈ 397093.1039551422
@test mean(P3) ≈ -243552.39197344126

#include("transient_creep.jl")
##@test mean(τ1) ≈ 136871.94383758004
#@test mean(τ2) ≈ 103128.05616241995




