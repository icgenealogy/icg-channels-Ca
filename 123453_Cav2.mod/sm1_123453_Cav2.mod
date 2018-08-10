NEURON
{
  SUFFIX Cav2 
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

  am = 0.18181808875052502     (/mV) 
  bm = -3.454533436400743     (1) 
  vhm = -53.99846304968979     (mV) 
  Am = 0.8675774266265222     (/ms) 
  b1m = 1.4478847724917445     (/mV) 
  c1m = -0.02024895838124459     (/mV2) 
  d1m = 6.971554967335881e-05     (/mV3) 
  b2m = -0.022473296756787352     (/mV) 
  c2m = 0.0014667715922738226     (/mV2) 
  d2m = -8.883067236822238e-06     (/mV3) 
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
  g = gbar*m*m
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