NEURON
{
  SUFFIX cat 
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

  ah = -0.09999995032102411     (/mV) 
  bh = 8.500030014266391     (1) 
  vhh = -78.0068284207718     (mV) 
  Ah = 63.02107505975269     (/ms) 
  b1h = 0.12321013264250948     (/mV) 
  c1h = 0.003325442452611093     (/mV2) 
  d1h = 7.105746535165396e-05     (/mV3) 
  b2h = 0.060408548390323134     (/mV) 
  c2h = -0.0005846468091863261     (/mV2) 
  d2h = 1.7336997231285762e-06     (/mV3) 

  am = 0.0999998522112253     (/mV) 
  bm = -5.9999909216144705     (1) 
  vhm = -40.48857616459581     (mV) 
  Am = 2.9467469423517136     (/ms) 
  b1m = -0.003828913345813563     (/mV) 
  c1m = -0.00039267373945284245     (/mV2) 
  d1m = 2.0197457395381623e-06     (/mV3) 
  b2m = -0.053860885850089996     (/mV) 
  c2m = -3.5293452709955836e-05     (/mV2) 
  d2m = 2.496397572229354e-06     (/mV3) 
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