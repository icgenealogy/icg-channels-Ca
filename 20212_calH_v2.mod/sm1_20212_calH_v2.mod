NEURON
{
  SUFFIX calH 
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

  am = 0.1505031914134968     (/mV) 
  bm = 0.8836067133181842     (1) 
  vhm = 6.4453818621427255     (mV) 
  Am = 14.113383941949275     (/ms) 
  b1m = -0.10717825170989581     (/mV) 
  c1m = 0.0003975772934714596     (/mV2) 
  d1m = 1.5107184505501452e-06     (/mV3) 
  b2m = -0.04512521319521392     (/mV) 
  c2m = 6.298332756666718e-05     (/mV2) 
  d2m = 4.305046359287878e-07     (/mV3) 
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