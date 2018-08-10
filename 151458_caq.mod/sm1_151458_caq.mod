NEURON
{
  SUFFIX caq 
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

  am = 0.15135218831198324     (/mV) 
  bm = -1.3621370644098325     (1) 
  vhm = -18.141537430666972     (mV) 
  Am = 0.7599999993241593     (/ms) 
  b1m = 1.65538767952024e-09     (/mV) 
  c1m = -1.1875502582163252e-08     (/mV2) 
  d1m = -7.105210228038413e-11     (/mV3) 
  b2m = 1.8060345780562148e-09     (/mV) 
  c2m = -1.1876807443458548e-08     (/mV2) 
  d2m = -7.109557587866617e-11     (/mV3) 
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