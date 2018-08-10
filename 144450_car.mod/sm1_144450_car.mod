NEURON
{
  SUFFIX car 
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

  ah = -1.0000946060570512     (/mV) 
  bh = 53.00548027645184     (1) 
  vhh = -55.679425747768015     (mV) 
  Ah = 8.003883873559953     (/ms) 
  b1h = 6.207478059387795e-05     (/mV) 
  c1h = -2.6291109586285995e-05     (/mV2) 
  d1h = 1.5880344245946716e-07     (/mV3) 
  b2h = 6.805620804628367e-05     (/mV) 
  c2h = -3.0169180410297357e-05     (/mV2) 
  d2h = 1.8238374932661076e-07     (/mV3) 

  am = 0.333213598906362     (/mV) 
  bm = -16.160688477197226     (1) 
  vhm = -14.797134498013657     (mV) 
  Am = 240.01580455838777     (/ms) 
  b1m = -0.013674755726736472     (/mV) 
  c1m = 9.430851745214936e-05     (/mV2) 
  d1m = -2.2731177400815642e-07     (/mV3) 
  b2m = -0.013659800466778435     (/mV) 
  c2m = -9.301096985956193e-05     (/mV2) 
  d2m = -2.1610150925707329e-07     (/mV3) 
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