NEURON
{
  SUFFIX Ca 
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

  ahCa = -0.11110705991131889     (/mV) 
  bhCa = -1.2222069670356654     (1) 
  vhhCa = 10.90123661720779     (mV) 
  AhCa = 100.00301503183957     (/ms) 
  b1hCa = -0.05521239723264206     (/mV) 
  c1hCa = -7.583060070220764e-06     (/mV2) 
  d1hCa = 6.316246399882316e-08     (/mV3) 
  b2hCa = -0.055810795394664534     (/mV) 
  c2hCa = -3.90922811522843e-06     (/mV2) 
  d2hCa = -9.956156880147007e-09     (/mV3) 

  amCa = 0.166669401075836     (/mV) 
  bmCa = -1.666742019931772     (1) 
  vhmCa = -3.8390123037318333     (mV) 
  AmCa = 9.44405861075215     (/ms) 
  b1mCa = -0.10228452380734493     (/mV) 
  c1mCa = -0.0025891012902730793     (/mV2) 
  d1mCa = 4.9359400704130744e-05     (/mV3) 
  b2mCa = -0.043545783253706394     (/mV) 
  c2mCa = 0.00036539656815636357     (/mV2) 
  d2mCa = 3.7513445374273977e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hCaInf 
  hCaTau 
  mCaInf 
  mCaTau 
}

STATE
{
  hCa
  mCa
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*hCa*mCa
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  hCa' = (hCaInf - hCa) / hCaTau 
  mCa' = (mCaInf - mCa) / mCaTau 
}

INITIAL
{
  rates(v)
  hCa = hCaInf 
  mCa = mCaInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hCaInf = 1/(1 + exp(-ahCa*v + bhCa)) 
    hCaTau = AhCa / ( exp(-(b1hCa*(v-vhhCa) + c1hCa*(v-vhhCa)^2 + d1hCa*(v-vhhCa)^3)) + exp((b2hCa*(v-vhhCa) + c2hCa*(v-vhhCa)^2 + d2hCa*(v-vhhCa)^3)) ) 

    mCaInf = 1/(1 + exp(-amCa*v + bmCa)) 
    mCaTau = AmCa / ( exp(-(b1mCa*(v-vhmCa) + c1mCa*(v-vhmCa)^2 + d1mCa*(v-vhmCa)^3)) + exp((b2mCa*(v-vhmCa) + c2mCa*(v-vhmCa)^2 + d2mCa*(v-vhmCa)^3)) ) 


  UNITSON
}