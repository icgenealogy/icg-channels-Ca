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

  ah = -1.989156754598786     (/mV) 
  bh = 81.56617420140077     (1) 
  vhh = -114.40779377120174     (mV) 
  Ah = 59.01453764138995     (/ms) 
  b1h = -0.0008969769114693284     (/mV) 
  c1h = 8.812634866735598e-06     (/mV2) 
  d1h = -2.3808589752549598e-08     (/mV3) 
  b2h = 0.0008937853250457099     (/mV) 
  c2h = -8.833807101961331e-06     (/mV2) 
  d2h = 2.3829550711832236e-08     (/mV3) 

  am = 0.9999986571776968     (/mV) 
  bm = -36.99990466483019     (1) 
  vhm = -133.76111083079616     (mV) 
  Am = 4.313581774964607     (/ms) 
  b1m = -0.003034885700176418     (/mV) 
  c1m = 1.6651946203941053e-05     (/mV2) 
  d1m = -2.9927598140897875e-08     (/mV3) 
  b2m = -0.11494957851013188     (/mV) 
  c2m = 0.0018716855287888216     (/mV2) 
  d2m = -1.2321652625858318e-05     (/mV3) 
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