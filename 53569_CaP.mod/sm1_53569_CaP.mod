NEURON
{
  SUFFIX CaP 
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

  am = 0.13513575480596637     (/mV) 
  bm = -7.027064281610257     (1) 
  vhm = -49.046128623290464     (mV) 
  Am = 17.065149426174994     (/ms) 
  b1m = -0.0021071550842488317     (/mV) 
  c1m = 1.631276799057114e-05     (/mV2) 
  d1m = 1.3359597891344154e-06     (/mV3) 
  b2m = 0.0020982572243393245     (/mV) 
  c2m = 8.939987624719191e-05     (/mV2) 
  d2m = -4.6182121287721886e-07     (/mV3) 
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