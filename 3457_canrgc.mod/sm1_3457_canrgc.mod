NEURON
{
  SUFFIX canrgc 
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

  ah = -0.15138513643948034     (/mV) 
  bh = 7.785881686527602     (1) 
  vhh = -41.842891226874485     (mV) 
  Ah = 86.09299570311524     (/ms) 
  b1h = 0.03732042839167298     (/mV) 
  c1h = -0.0009173017393380971     (/mV2) 
  d1h = 4.330080463900728e-06     (/mV3) 
  b2h = 0.07266140241586881     (/mV) 
  c2h = -0.0008866032332331805     (/mV2) 
  d2h = 3.2915512129686284e-06     (/mV3) 

  am = 0.12612287163170743     (/mV) 
  bm = -1.1382181574682104     (1) 
  vhm = -7.276298550321986     (mV) 
  Am = 5.897269988499524     (/ms) 
  b1m = 0.04896021607517551     (/mV) 
  c1m = -0.00011138853862021861     (/mV2) 
  d1m = -7.025070577442821e-07     (/mV3) 
  b2m = 0.07916223370360272     (/mV) 
  c2m = -0.0005961124872407353     (/mV2) 
  d2m = 1.7578372327347484e-06     (/mV3) 
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