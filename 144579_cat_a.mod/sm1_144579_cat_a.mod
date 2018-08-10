NEURON
{
  SUFFIX cat_a 
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

  ah = -0.1998972851378504     (/mV) 
  bh = 15.992365095772287     (1) 
  vhh = -73.29923113503305     (mV) 
  Ah = 269.060861154094     (/ms) 
  b1h = -0.15587001479795304     (/mV) 
  c1h = 0.0026854786212782266     (/mV2) 
  d1h = -1.072622946417937e-05     (/mV3) 
  b2h = -0.001542631469329915     (/mV) 
  c2h = 0.0004749288870404831     (/mV2) 
  d2h = -2.4371537818379723e-06     (/mV3) 

  am = 0.13512646557467972     (/mV) 
  bm = -7.026603160315603     (1) 
  vhm = -41.77857430412969     (mV) 
  Am = 5.003190129499428     (/ms) 
  b1m = -0.0012991174256840642     (/mV) 
  c1m = -0.000373485239134015     (/mV2) 
  d1m = 2.343820967085929e-06     (/mV3) 
  b2m = 0.0932851236055498     (/mV) 
  c2m = -0.001878838444131464     (/mV2) 
  d2m = 8.929820152175471e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}