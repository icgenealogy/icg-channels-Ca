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

  ah = -1.9914469255676504     (/mV) 
  bh = 81.6578667561202     (1) 
  vhh = -118.61769590876588     (mV) 
  Ah = 27.070746485084676     (/ms) 
  b1h = -0.005870433342897758     (/mV) 
  c1h = 3.679576154450658e-05     (/mV2) 
  d1h = -7.475023145497515e-08     (/mV3) 
  b2h = -0.11564432061707218     (/mV) 
  c2h = 0.002790020008088296     (/mV2) 
  d2h = -2.6076497955151125e-05     (/mV3) 

  am = 0.9999933444952398     (/mV) 
  bm = -37.6997221803186     (1) 
  vhm = -165.71287902255     (mV) 
  Am = 7.881822767428872     (/ms) 
  b1m = 4.953707453276648e-06     (/mV) 
  c1m = 0.0001029215658429196     (/mV2) 
  d1m = -2.1065843714545983e-07     (/mV3) 
  b2m = 0.008836393235002904     (/mV) 
  c2m = -3.383586457933087e-05     (/mV2) 
  d2m = 4.446775646576478e-08     (/mV3) 
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
  g = gbar*h*h*m*m*m
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