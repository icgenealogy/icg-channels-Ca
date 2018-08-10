NEURON
{
  SUFFIX it2 
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

  ah = -0.1999354221052158     (/mV) 
  bh = 15.99525020016996     (1) 
  vhh = -73.2866306763641     (mV) 
  Ah = 194.92260060878206     (/ms) 
  b1h = 0.00146294958777716     (/mV) 
  c1h = -0.0004753063200458634     (/mV2) 
  d1h = 2.4410151891540075e-06     (/mV3) 
  b2h = 0.15659101296886516     (/mV) 
  c2h = -0.002698915896550359     (/mV2) 
  d2h = 1.078192015340711e-05     (/mV3) 

  am = 0.1351348928927396     (/mV) 
  bm = -7.0270108791076575     (1) 
  vhm = -41.8393758572573     (mV) 
  Am = 1.8653154026644625     (/ms) 
  b1m = -0.0013005052187960598     (/mV) 
  c1m = -0.0003773228679674253     (/mV2) 
  d1m = 2.3653800474903694e-06     (/mV3) 
  b2m = 0.09400569557759016     (/mV) 
  c2m = -0.00188523365503835     (/mV2) 
  d2m = 8.94483933838267e-06     (/mV3) 
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