NEURON
{
  SUFFIX caL 
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

  am = 0.16393420920696708     (/mV) 
  bm = -5.737689827509348     (1) 
  vhm = -74.61998832140905     (mV) 
  Am = 0.22000000007231466     (/ms) 
  b1m = 6.236749821010841e-08     (/mV) 
  c1m = -3.523799547645604e-10     (/mV2) 
  d1m = 3.5082816464135165e-11     (/mV3) 
  b2m = 6.231729827335251e-08     (/mV) 
  c2m = -3.5017587322675286e-10     (/mV2) 
  d2m = 3.506587673046952e-11     (/mV3) 
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