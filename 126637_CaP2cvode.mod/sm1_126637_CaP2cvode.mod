NEURON
{
  SUFFIX CaP2cvode 
  USEION ca READ eca WRITE ica 
  RANGE gbar, g, ica
  GLOBAL eca
}

UNITS
{
  (S) = (siemens)
  (mV) = (millivolt)
  (mA) = (milliamp)
}

PARAMETER
{
  gbar = 1 (S/cm2)

  am = 0.13766654351676905     (/mV) 
  bm = -2.728868995749766     (1) 
  vhm = -19.366930599986556     (mV) 
  Am = 1.2279599177435863     (/ms) 
  b1m = 0.06882025678615114     (/mV) 
  c1m = 0.00014463186206630716     (/mV2) 
  d1m = -2.1897677450516944e-06     (/mV3) 
  b2m = 0.08269348383146904     (/mV) 
  c2m = -0.0009741551049549042     (/mV2) 
  d2m = 3.759375144173076e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}