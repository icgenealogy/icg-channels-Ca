NEURON
{
  SUFFIX ical3 
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

  am = 0.09962479281780093     (/mV) 
  bm = -1.9638110873714552     (1) 
  vhm = -33.69592105224357     (mV) 
  Am = 0.4077310335903505     (/ms) 
  b1m = -0.03615691375697529     (/mV) 
  c1m = 0.0003334823163726417     (/mV2) 
  d1m = -1.0237812392563076e-06     (/mV3) 
  b2m = -0.053927446926091756     (/mV) 
  c2m = -0.0007840914466289928     (/mV2) 
  d2m = -4.822756931623734e-06     (/mV3) 
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