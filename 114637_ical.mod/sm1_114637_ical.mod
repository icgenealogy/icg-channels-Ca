NEURON
{
  SUFFIX icalnew 
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

  am = 0.2500017006179861     (/mV) 
  bm = -6.2500374682971955     (1) 
  vhm = -93.84104742371152     (mV) 
  Am = 4.019999996405413     (/ms) 
  b1m = -1.6159877893627844e-10     (/mV) 
  c1m = -4.656949721755101e-10     (/mV2) 
  d1m = 3.780507896327917e-11     (/mV3) 
  b2m = -4.1766523575987943e-10     (/mV) 
  c2m = -4.6026126901511993e-10     (/mV2) 
  d2m = 3.7775132999452425e-11     (/mV3) 
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