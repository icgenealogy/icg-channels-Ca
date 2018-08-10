NEURON
{
  SUFFIX L_Ca 
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

  am = 0.15166231748011227     (/mV) 
  bm = -4.518223249507344     (1) 
  vhm = -630.8455685883267     (mV) 
  Am = 19.99619925027401     (/ms) 
  b1m = 0.9912075184570589     (/mV) 
  c1m = 0.006558265365234726     (/mV2) 
  d1m = -1.3328076120814432e-05     (/mV3) 
  b2m = -3.6768843112130477e-06     (/mV) 
  c2m = 6.526195527893806e-09     (/mV2) 
  d2m = -3.857836849675555e-12     (/mV3) 
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