NEURON
{
  SUFFIX Nca 
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

  ah = -0.046770649708366634     (/mV) 
  bh = 2.690880672362669     (1) 
  vhh = -7.51986962402407     (mV) 
  Ah = 153.17719792040774     (/ms) 
  b1h = -0.012631591191814966     (/mV) 
  c1h = 0.00020232411953144218     (/mV2) 
  d1h = 2.1288373169940657e-06     (/mV3) 
  b2h = 0.01263276797883468     (/mV) 
  c2h = 7.406006946140401e-05     (/mV2) 
  d2h = -8.644974868056153e-07     (/mV3) 

  am = 0.21890751547770643     (/mV) 
  bm = -7.232253903097635     (1) 
  vhm = -36.70502455911026     (mV) 
  Am = 4.309802123633813     (/ms) 
  b1m = 0.09486094588616899     (/mV) 
  c1m = 0.0012644288783173457     (/mV2) 
  d1m = 1.2993629597202561e-05     (/mV3) 
  b2m = 0.1367674850261803     (/mV) 
  c2m = -0.0016666869570061565     (/mV2) 
  d2m = 6.947830564249755e-06     (/mV3) 
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