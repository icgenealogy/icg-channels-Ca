NEURON
{
  SUFFIX caL 
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

  ah = -0.08400531763772376     (/mV) 
  bh = 1.1256679148799704     (1) 
  vhh = -48.42334587703105     (mV) 
  Ah = 29.54000002760079     (/ms) 
  b1h = 1.4944965025127885e-09     (/mV) 
  c1h = -7.874002394758541e-09     (/mV2) 
  d1h = -3.7446395204680715e-11     (/mV3) 
  b2h = 1.5809344030938598e-09     (/mV) 
  c2h = -7.872080706999396e-09     (/mV2) 
  d2h = -3.74875429686615e-11     (/mV3) 

  am = 0.14909792053238255     (/mV) 
  bm = -1.3269393148021411     (1) 
  vhm = 3.665308986583909     (mV) 
  Am = 0.5306445153853682     (/ms) 
  b1m = -0.058410889303705645     (/mV) 
  c1m = 0.0004640532770890054     (/mV2) 
  d1m = -2.3581308372780397e-06     (/mV3) 
  b2m = -0.02208870455985766     (/mV) 
  c2m = -2.4724271022453868e-05     (/mV2) 
  d2m = 3.0722322260087925e-07     (/mV3) 
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