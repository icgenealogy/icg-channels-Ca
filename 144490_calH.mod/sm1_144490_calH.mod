NEURON
{
  SUFFIX calH 
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

  ah = -1.9869072871243927     (/mV) 
  bh = 81.476185978848     (1) 
  vhh = -39.93349530282472     (mV) 
  Ah = 58.0201260278198     (/ms) 
  b1h = -0.008761566303428513     (/mV) 
  c1h = 3.903226698538039e-05     (/mV2) 
  d1h = -6.150475449786692e-08     (/mV3) 
  b2h = -0.008760149823132772     (/mV) 
  c2h = -3.7716850096959294e-05     (/mV2) 
  d2h = -5.111286580373653e-08     (/mV3) 

  am = 0.9999904286231086     (/mV) 
  bm = -36.99959762717022     (1) 
  vhm = -859.4743108628966     (mV) 
  Am = 0.08547103686483906     (/ms) 
  b1m = 0.7349695582048924     (/mV) 
  c1m = 0.0024127560734538843     (/mV2) 
  d1m = -3.846453471081479e-06     (/mV3) 
  b2m = -0.013851536771501102     (/mV) 
  c2m = 1.7056274251557422e-05     (/mV2) 
  d2m = -6.9882226760090666e-09     (/mV3) 
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
  g = gbar*h*m*m*m
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