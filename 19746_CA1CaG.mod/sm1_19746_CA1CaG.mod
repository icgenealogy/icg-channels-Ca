NEURON
{
  SUFFIX CA1CaG 
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

  ah = -0.13670106895873713     (/mV) 
  bh = 8.23850706613245     (1) 
  vhh = -61.215067571192314     (mV) 
  Ah = 14.139865669957793     (/ms) 
  b1h = -0.03290844180435014     (/mV) 
  c1h = -7.156690995835498e-05     (/mV2) 
  d1h = 6.622171501671847e-07     (/mV3) 
  b2h = -0.10509168148683236     (/mV) 
  c2h = -0.0001550570088793434     (/mV2) 
  d2h = -1.454618899398166e-06     (/mV3) 

  am = 0.14543351117347458     (/mV) 
  bm = -3.6591311854131696     (1) 
  vhm = -34.44416454783675     (mV) 
  Am = 1.6075529947917673     (/ms) 
  b1m = 0.0689209939017136     (/mV) 
  c1m = 0.0011934978542630932     (/mV2) 
  d1m = 8.325083273918954e-06     (/mV3) 
  b2m = 0.07917898134027455     (/mV) 
  c2m = -0.0007730562000175088     (/mV2) 
  d2m = 2.8360277076553285e-06     (/mV3) 
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