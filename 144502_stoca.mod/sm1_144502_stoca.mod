NEURON
{
  SUFFIX stoca 
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

  ah = -0.11624202409042367     (/mV) 
  bh = 9.939262268484514     (1) 
  vhh = -79.42947756671545     (mV) 
  Ah = 379.6214545987374     (/ms) 
  b1h = -0.09601644379135474     (/mV) 
  c1h = 0.0012834933424014558     (/mV2) 
  d1h = -4.5111022271343534e-06     (/mV3) 
  b2h = -0.009567315448149085     (/mV) 
  c2h = 0.0004773667208482246     (/mV2) 
  d2h = -2.155300662703629e-06     (/mV3) 

  am = 0.30300621955252993     (/mV) 
  bm = -16.682666710377244     (1) 
  vhm = -75.26966098881931     (mV) 
  Am = 164.84867571347226     (/ms) 
  b1m = -0.07435376992011108     (/mV) 
  c1m = 0.0013665505890111923     (/mV2) 
  d1m = -9.503073131300652e-06     (/mV3) 
  b2m = -1.21186446648282     (/mV) 
  c2m = 0.024862781406871207     (/mV2) 
  d2m = -1.3772212531838919e-05     (/mV3) 
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
  g = gbar*h*m
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