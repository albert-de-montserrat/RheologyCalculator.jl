= Tensor reduction operations

== Purely Maxwell Body

Model: visco-elastic

$bold(epsilon) = 1/(2 eta) bold(tau) + (bold(tau) - bold(tau)^o)/(2 G Delta t)$

First and final step:

$bold(epsilon) + (bold(tau)^o)/(2 G Delta t) = 1/(2 eta) bold(tau) + (bold(tau))/(2 G Delta t)$

Then:

$bold(epsilon)^("eff") = bold(epsilon) + (bold(tau)^o)/(2 G Delta t)$

== Kelvin-Voigt Body

Model: visco - (elastic-viscous).

Equations:

$bold(epsilon) = bold(epsilon)^(v_1) + bold(epsilon)^p = 1/(2 eta_1) bold(tau) + bold(epsilon)^p$

$bold(tau) = bold(tau)^(v_2) + bold(tau)^("el") = 2 eta_2 bold(epsilon)^p + (2 G Delta t bold(epsilon)^p + bold(tau^o))$

Rearranging the second equation gives:

$bold(tau) = (2 eta_2 + 2 G  Delta t) bold(epsilon)^p + bold(tau^o)$

$(bold(tau) - bold(tau^o)) 1 / (2 eta_2 + 2 G  Delta t) = (bold(tau) - bold(tau^o)) 1 / (2 eta^("eff")) = bold(epsilon)^p $

where 

$eta^("eff") = eta_2 + G Delta t$

In the most general case this would be:

$eta^("eff") = eta^("eff")_("KV") = sum^N_i eta_i$

where $eta^("eff")_("KV")$ is the effective viscosity of the Kelvin-Voigt element, $N$ is the number of branches in the Kelvin-Voigt element and $eta_i$ is the *effective* viscosity of each branch.

Plugging this into the first equation gives:

$bold(epsilon) = 1/(2 eta_1) bold(tau) + (bold(tau) - bold(tau^o)) 1 / (2 eta^("eff")_("KV"))$

and rearranging for an effective strain gives:

$bold(epsilon)^("eff") = bold(epsilon) + bold(tau^o) / (2 eta^("eff")_("KV")) = 1/(2 eta_1) bold(tau) + bold(tau) / (2 eta^("eff")_("KV"))$

== Mixed Kelvin-Voigt and Maxwell Body

Model: visco - (viscous-(elastic-viscous))

Equations:

$bold(epsilon) = bold(epsilon)^(v_1) + bold(epsilon)^p = 1/(2 eta_1) bold(tau) + bold(epsilon)^p$

$bold(tau) = bold(tau)^(v_2) + bold(tau)^("Max") = 2 eta_2 bold(epsilon)^p + bold(tau)^("Max")$

$bold(epsilon)^("Max") = bold(epsilon)^(p) =  bold(tau)^("Max") / (2 eta_3) +
(bold(tau)^("Max") - bold(tau)^(o))/ (2 G Delta t)$

Lets get $bold(tau)^("Max")$

$bold(epsilon)^p = bold(tau)^("Max") (1 / (2 eta_3) + 1/ (2 G Delta t)) -
 bold(tau)^(o)/ (2 G Delta t)$

$bold(tau)^("Max")  = (bold(epsilon)^p + bold(tau)^(o)/ (2 G Delta t)) (1 / (2 eta_3) + 1/ (2 G Delta t)) ^(-1) = 2 eta^("eff")_M (bold(epsilon)^p + bold(tau)^(o)/ (2 G Delta t)) $

where $eta^("eff")_M = 1 / (1 / (eta_3) + 1/ (G Delta t))$ is the Maxwell effective viscosity.

Plugging this into the second equation and solving for $bold(epsilon)^p$:

$bold(tau) = 2 eta_2 bold(epsilon)^p +  2 eta^("eff")_M (bold(epsilon)^p + bold(tau)^(o)/ (2 G Delta t))$

$bold(tau) = 2 eta_2 bold(epsilon)^p +  2 eta^("eff")_M bold(epsilon)^p + 2 eta^("eff")_M bold(tau)^(o)/ (2 G Delta t) $

$bold(epsilon)^p = (bold(tau) - 2 eta^("eff")_M bold(tau)^(o)/ (2 G Delta t))(2 eta_2 +  2 eta^("eff")_M)^(-1)$

Finally, plugging this into the first equation and rearranging for an effective strain gives:

$bold(epsilon) = 1/(2 eta_1) bold(tau) + (bold(tau) - 2 eta^("eff")_M bold(tau)^(o)/ (2 G Delta t))(2 eta_2 +  2 eta^("eff")_M)^(-1)$

= Generalized Maxwell Body

$
tau / (2 eta_("KV")) = epsilon^p + 1 / (2 eta_("KV")) (2 sum^N_i eta_(M_i)^star tau^o_i)
$

where $eta_("KV")$ is the effective viscosity of the Kelvin-Voigt element:

$
eta_("KV") = sum^N_i eta^("eff")_i
$

where $k$ is the branch number and $eta^("eff")_i$ is the effective viscosity of branch $i$; and $eta_(M_i)^star$ is the effective viscosity of the Maxwell element $i$:

$ 
eta_(M_i)^star = cases(
  eta_i / (eta_i + G_i Delta t) "if branch is Maxwell",
  1 "else",
) 
$

$
bold(r) = dot(bold(epsilon)) - f(bold(tau)) = 0
$ 
$
bold(tau)^(n+1) = bold(tau)^(n) - J^(-1) bold(r)
$ 

$cal(O)(4)$