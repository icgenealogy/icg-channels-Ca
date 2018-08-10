NEURON
{
  SUFFIX ical 
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

  ah = -0.04675414204348822     (/mV) 
  bh = 2.6906333907402864     (1) 
  vhh = -7.521430051013653     (mV) 
  Ah = 491.5884688368678     (/ms) 
  b1h = -0.012633724975852153     (/mV) 
  c1h = 0.00020229889312325483     (/mV2) 
  d1h = 2.128732242921264e-06     (/mV3) 
  b2h = 0.012633000677000345     (/mV) 
  c2h = 7.406851043817088e-05     (/mV2) 
  d2h = -8.645734596346878e-07     (/mV3) 

  am = 0.21890744451046024     (/mV) 
  bm = -7.232251286966226     (1) 
  vhm = -36.47581488603373     (mV) 
  Am = 13.735962339028566     (/ms) 
  b1m = 0.09285585760062998     (/mV) 
  c1m = 0.0011907123048920123     (/mV2) 
  d1m = 1.243172260648114e-05     (/mV3) 
  b2m = 0.13877083766620496     (/mV) 
  c2m = -0.0017094357174868467     (/mV2) 
  d2m = 7.2390569700017816e-06     (/mV3) 
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