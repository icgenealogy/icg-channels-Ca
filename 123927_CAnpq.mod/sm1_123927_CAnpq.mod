NEURON
{
  SUFFIX CAnpq 
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

  am = 0.12924973351316354     (/mV) 
  bm = -2.3249031522209993     (1) 
  vhm = 3.2218292385410456     (mV) 
  Am = 0.15360868731733404     (/ms) 
  b1m = 0.016650535709488924     (/mV) 
  c1m = -0.0002661373949476365     (/mV2) 
  d1m = 2.857270653203533e-08     (/mV3) 
  b2m = 0.11363079062367619     (/mV) 
  c2m = -0.0029340009607311205     (/mV2) 
  d2m = 4.329651840099537e-05     (/mV3) 
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