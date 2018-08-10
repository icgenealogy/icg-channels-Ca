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

  ah = -0.16149110250116724     (/mV) 
  bh = 10.835798884103498     (1) 
  vhh = -94.00250552542074     (mV) 
  Ah = 20.021414419885037     (/ms) 
  b1h = 0.0004338553416125248     (/mV) 
  c1h = -2.2905183215291864e-05     (/mV2) 
  d1h = 5.261995873945478e-08     (/mV3) 
  b2h = 0.00036582339387525916     (/mV) 
  c2h = -1.8369176619530805e-05     (/mV2) 
  d2h = -5.286242808886332e-09     (/mV3) 

  am = 0.1271377984389123     (/mV) 
  bm = -4.522308200045612     (1) 
  vhm = -16.406391267229846     (mV) 
  Am = 2.26826253029622     (/ms) 
  b1m = -0.007912655797785093     (/mV) 
  c1m = -0.00010950786697476281     (/mV2) 
  d1m = -2.3739383552960288e-08     (/mV3) 
  b2m = 0.0890021488919847     (/mV) 
  c2m = -0.0009181649540171597     (/mV2) 
  d2m = 1.4109094605125834e-06     (/mV3) 
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