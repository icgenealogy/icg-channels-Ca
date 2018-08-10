NEURON
{
  SUFFIX Ca_LVAst 
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

  ah = -0.15624681168442794     (/mV) 
  bh = 14.062261920210709     (1) 
  vhh = -31.443221162594888     (mV) 
  Ah = 16.215888697312508     (/ms) 
  b1h = -0.040269099633047624     (/mV) 
  c1h = 0.0016671285104727174     (/mV2) 
  d1h = -1.0188981019253005e-05     (/mV3) 
  b2h = 0.012002257884584939     (/mV) 
  c2h = 4.53382169364083e-05     (/mV2) 
  d2h = -7.825955184115957e-07     (/mV3) 

  am = 0.1666663734138382     (/mV) 
  bm = -6.666649139870156     (1) 
  vhm = -25.263726294314488     (mV) 
  Am = 5.925889240735507     (/ms) 
  b1m = -0.10102327844074367     (/mV) 
  c1m = 0.0027436134461012074     (/mV2) 
  d1m = -1.5018380882248227e-05     (/mV3) 
  b2m = 0.013555740626518788     (/mV) 
  c2m = 6.98571614374866e-05     (/mV2) 
  d2m = -8.026671851405405e-07     (/mV3) 
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