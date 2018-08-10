NEURON
{
  SUFFIX cal 
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

  am = 0.12932456725317815     (/mV) 
  bm = -2.5764776391463147     (1) 
  vhm = -15.483311867141992     (mV) 
  Am = 4.2004657745672604     (/ms) 
  b1m = 0.05538211280620096     (/mV) 
  c1m = 0.000653512171619316     (/mV2) 
  d1m = 3.3903195574316375e-06     (/mV3) 
  b2m = 0.07366252724619504     (/mV) 
  c2m = -0.0009477441611305237     (/mV2) 
  d2m = 3.898954399011863e-06     (/mV3) 
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