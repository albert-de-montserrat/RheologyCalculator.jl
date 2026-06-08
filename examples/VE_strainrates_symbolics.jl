using Symbolics

@variables εxx,εxz, G, τxx, τxz, η, τxx_o, τxz_o, Δt


# Serial case:
# --⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--
r = [εxx, εxz] - 1/2η*[τxx,τxz] - [(τxx-τxx_o),(τxz-τxz_o)] *1/(2G*Δt)
# εxx_eff = εxx + τxx_o/(2G*Δt)
# εxz_eff = εxz + τxz_o/(2G*Δt)
# εII_eff = (1/2η+ 1/(2G*Δt)) τII


# Parallel case:
# |--⟦▪̲̅▫̲̅▫̲̅▫̲̅¹--|
# |--/\/\/¹--|
r_ε =  [εxx, εxz] 
r_τ = -[τxx,τxz] + (2η + 2G*Δt)*[εxx, εxz]  +  [τxx_o,τxz_o]
# rII = -τII + (2η + 2G*Δt)*εII + τII_o

# 2 parallel blocks in series:
# |--⟦▪̲̅▫̲̅▫̲̅▫̲̅--|     |--⟦▪̲̅▫̲̅▫̲̅▫̲̅--|
# |--/\/\/--|  +  |--/\/\/--| 
r   = r_ε - 1/2η*r_τ - [(τxx-τxx_o),(τxz-τxz_o)] *1/(2G*Δt)
r_ε =   [εxx, εxz] - [εxx1, εxz1] - [εxx2, εxz2]
r_τ =  -[τxx,τxz] + (2η1 + 2G1*Δt)*[εxx1, εxz1]  +  [τxx1_o,τxz1_o] + (2η2 + 2G2*Δt)*[εxx2, εxz2]  +  [τxx2_o,τxz2_o]

#r_εII =  εII - εII1 - εII2
#r_τ   =  -τII + (2η1 + 2G1*Δt)*εII1  +  τII_eff_o + (2η2 + 2G2*Δt)*[εxx2, εxz2]  
# with:
#   τxx_eff_o = τxx1_o + τxx2_o
#   τxz_eff_o = τxz1_o + τxz2_o
#   τII_eff_o = sqrt(τxx_eff_o^2 + τxz_eff_o^2)

