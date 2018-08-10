NEURON
{
  SUFFIX cal_gp 
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

  am = 0.08333339554523847     (/mV) 
  bm = -3.3333369973209708     (1) 
  vhm = -1.9970074898590324     (mV) 
  Am = 15.601309899724844     (/ms) 
  b1m = -0.08305724845192641     (/mV) 
  c1m = -2.318119386869201e-05     (/mV2) 
  d1m = 5.457395181377827e-07     (/mV3) 
  b2m = -0.09055220658286599     (/mV) 
  c2m = 2.504142472933387e-05     (/mV2) 
  d2m = 6.460480503570698e-07     (/mV3) 
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