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

  am = 0.20124458401935436     (/mV) 
  bm = -5.977196549255123     (1) 
  vhm = -35.41611417840934     (mV) 
  Am = 1.6579801640524516     (/ms) 
  b1m = 0.10567645077326601     (/mV) 
  c1m = 0.001594186128887317     (/mV2) 
  d1m = 1.568131567759093e-05     (/mV3) 
  b2m = 0.10883673182322906     (/mV) 
  c2m = -0.0012217133930393195     (/mV2) 
  d2m = 4.823520182060073e-06     (/mV3) 
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