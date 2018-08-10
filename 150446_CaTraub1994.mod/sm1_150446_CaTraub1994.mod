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

  am = 0.1293113692464191     (/mV) 
  bm = -3.223095409425029     (1) 
  vhm = -22.631709985642996     (mV) 
  Am = 4.2586405779485625     (/ms) 
  b1m = 0.06152991613567684     (/mV) 
  c1m = 0.0008327725968501517     (/mV2) 
  d1m = 4.842896788415979e-06     (/mV3) 
  b2m = 0.0692511737818668     (/mV) 
  c2m = -0.000820218976642595     (/mV2) 
  d2m = 3.1281820296435065e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}