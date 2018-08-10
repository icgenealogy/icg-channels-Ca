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

  ah = -1.9891567068720388     (/mV) 
  bh = 81.56617229233139     (1) 
  vhh = -114.40779377120174     (mV) 
  Ah = 59.01453764138995     (/ms) 
  b1h = -0.0008969769114693284     (/mV) 
  c1h = 8.812634866735598e-06     (/mV2) 
  d1h = -2.3808589752549598e-08     (/mV3) 
  b2h = 0.0008937853250457099     (/mV) 
  c2h = -8.833807101961331e-06     (/mV2) 
  d2h = 2.3829550711832236e-08     (/mV3) 

  am = 0.9999986571776309     (/mV) 
  bm = -36.99990466482783     (1) 
  vhm = -133.56930613375994     (mV) 
  Am = 4.144679461820126     (/ms) 
  b1m = 0.13242172418260698     (/mV) 
  c1m = -0.0022044614125701226     (/mV2) 
  d1m = 1.4350823248911082e-05     (/mV3) 
  b2m = 0.0023803206102484607     (/mV) 
  c2m = -1.3350032953403418e-05     (/mV2) 
  d2m = 2.443723026999785e-08     (/mV3) 
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