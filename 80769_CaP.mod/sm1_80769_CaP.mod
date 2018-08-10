NEURON
{
  SUFFIX CaP 
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

  am = 0.18181801903905817     (/mV) 
  bm = -3.4545356890722987     (1) 
  vhm = -52.914796143261974     (mV) 
  Am = 0.742453497302555     (/ms) 
  b1m = 0.566606407681698     (/mV) 
  c1m = 0.03243010407461884     (/mV2) 
  d1m = 0.000476049464435519     (/mV3) 
  b2m = -0.018851548125464156     (/mV) 
  c2m = 0.0013270456353807326     (/mV2) 
  d2m = -7.1876308393655745e-06     (/mV3) 
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