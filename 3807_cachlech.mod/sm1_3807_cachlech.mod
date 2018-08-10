NEURON
{
  SUFFIX cach 
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

  am = 0.26152461881288475     (/mV) 
  bm = -1.5207487056499362     (1) 
  vhm = -4.859649441923712     (mV) 
  Am = 0.6782780447532111     (/ms) 
  b1m = 0.09261194591707365     (/mV) 
  c1m = 0.00018083258321583906     (/mV2) 
  d1m = 3.864431183185858e-06     (/mV3) 
  b2m = 0.1912906355999508     (/mV) 
  c2m = -0.0036323577086408688     (/mV2) 
  d2m = 4.870871446077021e-05     (/mV3) 
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