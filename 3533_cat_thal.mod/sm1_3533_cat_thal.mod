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

  ah = -0.0847455686712694     (/mV) 
  bh = 5.974567559017242     (1) 
  vhh = -64.98859200020733     (mV) 
  Ah = 19.251939809132747     (/ms) 
  b1h = -0.04975393862791725     (/mV) 
  c1h = -1.0112447949273394e-05     (/mV2) 
  d1h = 1.1958825387375204e-07     (/mV3) 
  b2h = -0.049701167445757075     (/mV) 
  c2h = 1.127200743738595e-05     (/mV2) 
  d2h = 1.4070439525222663e-07     (/mV3) 

  am = 0.07352925833823976     (/mV) 
  bm = -1.6249922485766928     (1) 
  vhm = -18.787095696261215     (mV) 
  Am = 1.9381933056309657     (/ms) 
  b1m = -0.05377312303261206     (/mV) 
  c1m = 7.778235162552661e-05     (/mV2) 
  d1m = -2.668941308048556e-07     (/mV3) 
  b2m = -0.04773965198276538     (/mV) 
  c2m = 2.9038095104184296e-05     (/mV2) 
  d2m = 3.263327566945867e-07     (/mV3) 
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
  g = gbar*h*m
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