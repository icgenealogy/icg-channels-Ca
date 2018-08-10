NEURON
{
  SUFFIX captain 
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

  am = 0.1818180537007923     (/mV) 
  bm = -3.4545349952189013     (1) 
  vhm = -52.96424003255698     (mV) 
  Am = 3.7030082917787177     (/ms) 
  b1m = 0.023177767480888154     (/mV) 
  c1m = -0.001445621503734246     (/mV2) 
  d1m = 7.776796605888221e-06     (/mV3) 
  b2m = -0.6077466883094232     (/mV) 
  c2m = -0.03721396354953574     (/mV2) 
  d2m = -0.0005849708825679597     (/mV3) 
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