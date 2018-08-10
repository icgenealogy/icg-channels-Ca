NEURON
{
  SUFFIX it2 
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

  ah = -0.19993641447300078     (/mV) 
  bh = 15.59545185278995     (1) 
  vhh = -84.69319105929424     (mV) 
  Ah = 324.40221285121385     (/ms) 
  b1h = -0.1072236806006989     (/mV) 
  c1h = 0.0011107830626456536     (/mV2) 
  d1h = -3.4005490961863287e-06     (/mV3) 
  b2h = -0.14580630240655806     (/mV) 
  c2h = -0.011270563014460005     (/mV2) 
  d2h = -0.0003600332170470643     (/mV3) 

  am = 0.1351349885896036     (/mV) 
  bm = -6.756746193023653     (1) 
  vhm = -39.28513739369515     (mV) 
  Am = 1.831363141931621     (/ms) 
  b1m = -0.0019295727910947664     (/mV) 
  c1m = -0.0003586802967310402     (/mV2) 
  d1m = 2.2887167444194126e-06     (/mV3) 
  b2m = 0.09297034311990654     (/mV) 
  c2m = -0.0019036850060264571     (/mV2) 
  d2m = 9.167748757423777e-06     (/mV3) 
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