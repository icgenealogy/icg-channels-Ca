NEURON
{
  SUFFIX iCaL 
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

  am = 0.10000007170344094     (/mV) 
  bm = -0.9999971095854324     (1) 
  vhm = -7.078830964792193     (mV) 
  Am = 0.6691793204935811     (/ms) 
  b1m = -0.06207755734408183     (/mV) 
  c1m = 0.0008336720687311614     (/mV2) 
  d1m = -3.6175279095092527e-06     (/mV3) 
  b2m = -0.03641270064457924     (/mV) 
  c2m = -0.00028131284036047417     (/mV2) 
  d2m = -9.783769916554986e-07     (/mV3) 
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