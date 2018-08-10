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

  ah = -0.24991563192038063     (/mV) 
  bh = 19.993520461103948     (1) 
  vhh = -86.89423260012055     (mV) 
  Ah = 198.50134927840008     (/ms) 
  b1h = 0.15962630786019671     (/mV) 
  c1h = 0.015537876297460255     (/mV2) 
  d1h = 0.0006132962349555777     (/mV3) 
  b2h = 0.10113962310867243     (/mV) 
  c2h = -0.0009703746263678161     (/mV2) 
  d2h = 2.8299167720412696e-06     (/mV3) 

  am = 0.1612768136846537     (/mV) 
  bm = -9.031522495665614     (1) 
  vhm = -65.54585678072567     (mV) 
  Am = 8.296427098763793     (/ms) 
  b1m = -0.07649122297760175     (/mV) 
  c1m = 0.0005184883115865568     (/mV2) 
  d1m = -1.1307597229051107e-06     (/mV3) 
  b2m = -0.023972230766038577     (/mV) 
  c2m = 0.00035361196857312     (/mV2) 
  d2m = -1.8431648546493287e-06     (/mV3) 
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