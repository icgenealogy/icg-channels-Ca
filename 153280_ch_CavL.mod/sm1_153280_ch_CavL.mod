NEURON
{
  SUFFIX ch_CavL 
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

  am = 0.18003452589837368     (/mV) 
  bm = -0.2399761101046787     (1) 
  vhm = -1.1426119281236924     (mV) 
  Am = 3.0585886221902183     (/ms) 
  b1m = 0.09080800139402823     (/mV) 
  c1m = -3.622565751396066e-05     (/mV2) 
  d1m = -1.283079439292832e-06     (/mV3) 
  b2m = 0.0896672552529819     (/mV) 
  c2m = -0.0001354760059778664     (/mV2) 
  d2m = -9.601185612945335e-07     (/mV3) 
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