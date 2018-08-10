NEURON
{
  SUFFIX cal 
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

  am = 0.17994971830878181     (/mV) 
  bm = -0.2389133252574028     (1) 
  vhm = -90.87089418072004     (mV) 
  Am = 0.7305294298715352     (/ms) 
  b1m = 0.0565831204631381     (/mV) 
  c1m = 0.0004417160745277329     (/mV2) 
  d1m = -4.118561738668882e-06     (/mV3) 
  b2m = 0.0022312739538111943     (/mV) 
  c2m = 0.0003439240456879204     (/mV2) 
  d2m = -1.4391654276682166e-06     (/mV3) 
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