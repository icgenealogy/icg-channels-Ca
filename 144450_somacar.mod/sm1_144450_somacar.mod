NEURON
{
  SUFFIX somacar 
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

  ah = -1.0006913839724847     (/mV) 
  bh = 62.046185013965705     (1) 
  vhh = -126.05283695746009     (mV) 
  Ah = 5.02906691835457     (/ms) 
  b1h = -7.239789706064386e-05     (/mV) 
  c1h = 4.448727375886119e-07     (/mV2) 
  d1h = -8.848853240677416e-10     (/mV3) 
  b2h = -0.4351462455999367     (/mV) 
  c2h = 0.009316420041156941     (/mV2) 
  d2h = -6.076288122405101e-05     (/mV3) 

  am = 0.333216196998508     (/mV) 
  bm = -19.993806508545866     (1) 
  vhm = -18.865060693366782     (mV) 
  Am = 200.00755905746686     (/ms) 
  b1m = -0.013993303120205598     (/mV) 
  c1m = 9.8517649722763e-05     (/mV2) 
  d1m = -2.415866923072075e-07     (/mV3) 
  b2m = -0.01397401027785781     (/mV) 
  c2m = -9.788290866607856e-05     (/mV2) 
  d2m = -2.3807148884426978e-07     (/mV3) 
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