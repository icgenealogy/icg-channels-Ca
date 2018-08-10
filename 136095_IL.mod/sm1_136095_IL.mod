NEURON
{
  SUFFIX ical 
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

  am = 0.09962479312235104     (/mV) 
  bm = -1.9638111021163918     (1) 
  vhm = -26.21420811932512     (mV) 
  Am = 0.41392309299550895     (/ms) 
  b1m = -0.04510500948337267     (/mV) 
  c1m = 0.0004960967857861752     (/mV2) 
  d1m = -1.767638066153981e-06     (/mV3) 
  b2m = -0.05180689684197032     (/mV) 
  c2m = -0.0007737499169737841     (/mV2) 
  d2m = -5.097573168935704e-06     (/mV3) 
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
  g = gbar*m*m
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