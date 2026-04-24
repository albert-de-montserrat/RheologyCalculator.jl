#set math.equation(numbering: "(1)") 

#show ref: it => {
  let eq = math.equation
  let el = it.element
  // Skip all other references.
  if el == none or el.func() != eq { return it }
  // Override equation references.
  link(el.location(), numbering(
    el.numbering,
    ..counter(eq).at(el.location())
  ))
}

= Generalized Maxwell Body

$
  epsilon = epsilon^p
$<eq1>
$
  tau^p = tau = tau_1 + tau_2 + ... + tau_n = sum^("NK")_i tau^p_i
$<eq2>
$
  epsilon_i = sum^("NM"_i)_j epsilon_i_j
$<eq3>

where NK is the number of Kelvin elements, and $"NM"_i$ is the number of Maxwell elements in the i-th Kelvin element.

For the general visco-elastic Maxwell model we then have

$
  epsilon_i = (1) / (2 eta_i) tau_i + 1 / (2 G_i Delta t) (tau_i - tau_i^o)
$

, and solving for $tau_i$ gives

$
  tau_i = (epsilon_i + tau_i^o / (2 G_i Delta t)) (1 / (2 eta_i) + 1 / (2 G_i Delta t))^(-1) = (epsilon_i + tau_i^o / (2 G_i Delta t)) 2 eta_M_i
$<eq5>

where $eta_M_i$ is the Maxwell viscosity of the i-th Kelvin element, defined as
$
  eta_M_i = 1 / (eta_i) + 1 / (G_i Delta t)
$

Now, plugging eq. @eq5 in eq. @eq2 gives
$
  tau = (epsilon_1 + tau_1^o / (2 G_1 Delta t)) 2 eta_M_1 + (epsilon_2 + tau_2^o / (2 G_2 Delta t)) 2 eta_M_2 + ... + (epsilon_n + tau_n^o / (2 G_n Delta t)) 2 eta_M_n
$

$
  tau = 2 eta_K epsilon + sum^("NK")_i eta_M_i tau_i^o / (G_i Delta t)
$
$
  tau / (2 eta_K) =  epsilon + (1)/(2 eta_K) sum^("NK")_i eta_M_i tau_i^o / (G_i Delta t)
$
where $eta_K = sum^("NK")_i eta_M_i$ is the effective viscosity of the Kelvin element.
