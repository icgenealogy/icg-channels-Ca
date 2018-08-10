NEURON
{
  SUFFIX iT 
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

  ah = -0.24810740998630604     (/mV) 
  bh = 20.841216437958977     (1) 
  vhh = -88.23772823751112     (mV) 
  Ah = 147.50360615158505     (/ms) 
  b1h = 0.14915453218925365     (/mV) 
  c1h = 0.014323724470285707     (/mV2) 
  d1h = 0.0006031407331160682     (/mV3) 
  b2h = 0.09902767351645209     (/mV) 
  c2h = -0.000936495968917972     (/mV2) 
  d2h = 2.6975795444527994e-06     (/mV3) 

  am = 0.16128987149069665     (/mV) 
  bm = -9.758036622450902     (1) 
  vhm = -67.21293874917212     (mV) 
  Am = 4.560199585095358     (/ms) 
  b1m = -0.0743466723641632     (/mV) 
  c1m = 0.0004645622301415372     (/mV2) 
  d1m = -7.792363214331785e-07     (/mV3) 
  b2m = -0.026115133463849194     (/mV) 
  c2m = 0.0003843140889486081     (/mV2) 
  d2m = -7.510467617320505e-07     (/mV3) 
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