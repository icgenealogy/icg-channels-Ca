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

  ah = -0.1993701102243082     (/mV) 
  bh = 15.552607088400995     (1) 
  vhh = -70.74797978110308     (mV) 
  Ah = 780.0202214822399     (/ms) 
  b1h = -0.0007239347893162163     (/mV) 
  c1h = -0.0004492050253954531     (/mV2) 
  d1h = 2.389316210279974e-06     (/mV3) 
  b2h = 0.1591755851280633     (/mV) 
  c2h = -0.0028332638975581605     (/mV2) 
  d2h = 1.1557662335472844e-05     (/mV3) 

  am = 0.13503088765279433     (/mV) 
  bm = -6.751568694847954     (1) 
  vhm = -27.26821236827414     (mV) 
  Am = 8.183313561098226     (/ms) 
  b1m = -0.024491769701262347     (/mV) 
  c1m = 0.0007401748714216147     (/mV2) 
  d1m = -4.369811609741751e-06     (/mV3) 
  b2m = 0.006225975713938677     (/mV) 
  c2m = 0.0001441075815652667     (/mV2) 
  d2m = -1.2821604127986614e-06     (/mV3) 
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