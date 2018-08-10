NEURON
{
  SUFFIX cal_ch 
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

  am = 0.17857145362324947     (/mV) 
  bm = -0.2678423070687405     (1) 
  vhm = -2.574835522165784     (mV) 
  Am = 3.000845427899258     (/ms) 
  b1m = -0.07876840829765058     (/mV) 
  c1m = -0.00012153801250790632     (/mV2) 
  d1m = 1.3696029495518017e-06     (/mV3) 
  b2m = -0.09474670635626732     (/mV) 
  c2m = -0.00010674301055119421     (/mV2) 
  d2m = -1.7498744273957993e-07     (/mV3) 
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