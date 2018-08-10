NEURON
{
  SUFFIX cap_gp 
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

  am = 0.1864634979454993     (/mV) 
  bm = -3.643733100587134     (1) 
  vhm = -53.22854318861984     (mV) 
  Am = 20.01371152136168     (/ms) 
  b1m = 0.5550975346015757     (/mV) 
  c1m = 0.03347031835477279     (/mV2) 
  d1m = 0.0005314164706742988     (/mV3) 
  b2m = -0.019602327428180837     (/mV) 
  c2m = 0.0013804558993759492     (/mV2) 
  d2m = -7.474463922559169e-06     (/mV3) 
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