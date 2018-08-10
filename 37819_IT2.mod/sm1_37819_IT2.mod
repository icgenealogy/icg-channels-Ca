NEURON
{
  SUFFIX itre 
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

  ah = -0.19991883832503526     (/mV) 
  bh = 15.994033779514448     (1) 
  vhh = -73.28721768064244     (mV) 
  Ah = 247.08041609696795     (/ms) 
  b1h = -0.15656979726963227     (/mV) 
  c1h = 0.002698379410248853     (/mV2) 
  d1h = -1.0779545661705485e-05     (/mV3) 
  b2h = -0.0014653626460391673     (/mV) 
  c2h = 0.00047535381948576037     (/mV2) 
  d2h = -2.4412032775394574e-06     (/mV3) 

  am = 0.1351348928927435     (/mV) 
  bm = -7.027010879107868     (1) 
  vhm = -41.74532953569871     (mV) 
  Am = 1.526792489608477     (/ms) 
  b1m = -0.0012933630045428642     (/mV) 
  c1m = -0.00037246463054311     (/mV2) 
  d1m = 2.3392969049301247e-06     (/mV3) 
  b2m = 0.09259337907813635     (/mV) 
  c2m = -0.0018670898787190508     (/mV2) 
  d2m = 8.881074003635688e-06     (/mV3) 
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