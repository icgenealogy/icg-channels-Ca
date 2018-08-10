NEURON
{
  SUFFIX sca 
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

  ah = -0.04675458451872496     (/mV) 
  bh = 2.69058230908676     (1) 
  vhh = -7.5103599235340335     (mV) 
  Ah = 459.4626458923859     (/ms) 
  b1h = -0.012631440037052204     (/mV) 
  c1h = 0.0002023686947209345     (/mV2) 
  d1h = 2.128195943084212e-06     (/mV3) 
  b2h = 0.012632354795081607     (/mV) 
  c2h = 7.404574993297515e-05     (/mV2) 
  d2h = -8.644378457563955e-07     (/mV3) 

  am = 0.21890748760109224     (/mV) 
  bm = -7.232252874417753     (1) 
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