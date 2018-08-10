NEURON
{
  SUFFIX ca 
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

  ah = -0.0467658618593729     (/mV) 
  bh = 2.690599598234306     (1) 
  vhh = -7.518798921183556     (mV) 
  Ah = 153.18308397675926     (/ms) 
  b1h = -0.012631073440609252     (/mV) 
  c1h = 0.00020222775285743134     (/mV2) 
  d1h = 2.127978505344759e-06     (/mV3) 
  b2h = 0.012631058218590683     (/mV) 
  c2h = 7.405818944391527e-05     (/mV2) 
  d2h = -8.643063740142366e-07     (/mV3) 

  am = 0.2184341731302147     (/mV) 
  bm = -7.216316266832218     (1) 
  vhm = -36.60782636227561     (mV) 
  Am = 4.292185112542166     (/ms) 
  b1m = -0.13703238153257966     (/mV) 
  c1m = 0.0016745436991185958     (/mV2) 
  d1m = -6.995901086942127e-06     (/mV3) 
  b2m = -0.09360691599088436     (/mV) 
  c2m = -0.001218310871591697     (/mV2) 
  d2m = -1.2495245948411466e-05     (/mV3) 
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