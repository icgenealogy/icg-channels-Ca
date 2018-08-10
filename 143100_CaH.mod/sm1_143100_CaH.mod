NEURON
{
  SUFFIX CaH 
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

  am = 0.1428578451552806     (/mV) 
  bm = -2.8571637362318403     (1) 
  vhm = -92.10267662267806     (mV) 
  Am = 0.41999999979797964     (/ms) 
  b1m = -1.4649774033364614e-09     (/mV) 
  c1m = 1.3809326284340482e-10     (/mV2) 
  d1m = -3.4489615888870275e-11     (/mV3) 
  b2m = -1.6547046011454533e-09     (/mV) 
  c2m = 1.4244725775976686e-10     (/mV2) 
  d2m = -3.451477716898224e-11     (/mV3) 
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