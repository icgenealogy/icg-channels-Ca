NEURON
{
  SUFFIX CaP 
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

  ah = -0.1403232246596755     (/mV) 
  bh = 4.924035291552498     (1) 
  vhh = -50.48157246001848     (mV) 
  Ah = 1291.8603233936963     (/ms) 
  b1h = 0.053026529107413095     (/mV) 
  c1h = 0.0014781258733998872     (/mV2) 
  d1h = 1.377600766463382e-05     (/mV3) 
  b2h = 0.06787110416301069     (/mV) 
  c2h = -0.000722574399027309     (/mV2) 
  d2h = 2.4099172797485113e-06     (/mV3) 

  am = 0.13766289293280304     (/mV) 
  bm = -2.728731202457768     (1) 
  vhm = -19.366930599986556     (mV) 
  Am = 1.2279599177435863     (/ms) 
  b1m = 0.06882025678615114     (/mV) 
  c1m = 0.00014463186206630716     (/mV2) 
  d1m = -2.1897677450516944e-06     (/mV3) 
  b2m = 0.08269348383146904     (/mV) 
  c2m = -0.0009741551049549042     (/mV2) 
  d2m = 3.759375144173076e-06     (/mV3) 
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