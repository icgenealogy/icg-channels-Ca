NEURON
{
  SUFFIX Ca_LVA 
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

  ah = -0.15624763504024078     (/mV) 
  bh = 14.062321598361498     (1) 
  vhh = -31.3999074293399     (mV) 
  Ah = 12.621162689720823     (/ms) 
  b1h = -0.012008470424073693     (/mV) 
  c1h = -4.5183955287229824e-05     (/mV2) 
  d1h = 7.823566161688848e-07     (/mV3) 
  b2h = 0.04011542272068276     (/mV) 
  c2h = -0.0016652769486443526     (/mV2) 
  d2h = 1.01860181437952e-05     (/mV3) 

  am = 0.16666637341392904     (/mV) 
  bm = -6.666649139873668     (1) 
  vhm = -50.32482503206272     (mV) 
  Am = 12.977013405935857     (/ms) 
  b1m = 0.06394290240394067     (/mV) 
  c1m = 0.001966999271368943     (/mV2) 
  d1m = 1.9383617731364478e-05     (/mV3) 
  b2m = 0.09066209354543144     (/mV) 
  c2m = -0.001049840467222675     (/mV2) 
  d2m = 3.7098834170730156e-06     (/mV3) 
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