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

  ah = -0.13909052281319922     (/mV) 
  bh = 11.893874964413028     (1) 
  vhh = -91.7620439065106     (mV) 
  Ah = 1588.1488505099746     (/ms) 
  b1h = 0.140855640183455     (/mV) 
  c1h = 0.01926002159351353     (/mV2) 
  d1h = 0.0013766419210977207     (/mV3) 
  b2h = 0.18479217438471815     (/mV) 
  c2h = -0.0018751144065440985     (/mV2) 
  d2h = 5.608442757528087e-06     (/mV3) 

  am = 0.13280205068927886     (/mV) 
  bm = -7.883128185126603     (1) 
  vhm = -76.15794193958664     (mV) 
  Am = 2.8751012238774774     (/ms) 
  b1m = -0.07848183934066535     (/mV) 
  c1m = 0.0007385636271019548     (/mV2) 
  d1m = -2.158742533737386e-06     (/mV3) 
  b2m = -0.08881754997164998     (/mV) 
  c2m = -0.003083138628650478     (/mV2) 
  d2m = -5.7053540387437825e-05     (/mV3) 
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