NEURON
{
  SUFFIX cat 
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

  ah = -0.2702041800675826     (/mV) 
  bh = 19.45493031801708     (1) 
  vhh = -76.85012172147097     (mV) 
  Ah = 191.47515882901646     (/ms) 
  b1h = 0.08971510587307359     (/mV) 
  c1h = 0.0051454151187288865     (/mV2) 
  d1h = 9.773084021043759e-05     (/mV3) 
  b2h = 0.1813315715387304     (/mV) 
  c2h = -0.0034169617810854013     (/mV2) 
  d2h = 2.4875861971973604e-05     (/mV3) 

  am = 0.0838229282668678     (/mV) 
  bm = -3.5359316015919857     (1) 
  vhm = -30.18320887873785     (mV) 
  Am = 1.4001104469526071     (/ms) 
  b1m = -0.001879432694929305     (/mV) 
  c1m = -0.0002554044182165934     (/mV2) 
  d1m = 1.8615738070051712e-06     (/mV3) 
  b2m = 0.13703153408751398     (/mV) 
  c2m = -0.003348413425202571     (/mV2) 
  d2m = 1.7484501691508775e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}