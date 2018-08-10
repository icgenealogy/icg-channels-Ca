NEURON
{
  SUFFIX car 
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

  ah = -1.0000107239303506     (/mV) 
  bh = 50.00067610502119     (1) 
  vhh = -111.38191048516629     (mV) 
  Ah = 40.03728823246728     (/ms) 
  b1h = -3.2843523895873505e-05     (/mV) 
  c1h = 3.541030067158391e-07     (/mV2) 
  d1h = -1.010900774864131e-09     (/mV3) 
  b2h = 3.154019382892582e-05     (/mV) 
  c2h = -3.4714398821208883e-07     (/mV2) 
  d2h = 1.0034419385509805e-09     (/mV3) 

  am = 0.33326109160632683     (/mV) 
  bm = -14.496800321715583     (1) 
  vhm = -14.21063984217272     (mV) 
  Am = 139.99747303132634     (/ms) 
  b1m = -0.014158801310273758     (/mV) 
  c1m = 0.00010096633940312131     (/mV2) 
  d1m = -2.510913041835608e-07     (/mV3) 
  b2m = -0.014148818332125354     (/mV) 
  c2m = -0.00010010378003857715     (/mV2) 
  d2m = -2.43458586985669e-07     (/mV3) 
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