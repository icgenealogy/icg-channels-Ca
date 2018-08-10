NEURON
{
  SUFFIX CaLVA 
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

  ah = -0.24997086297611462     (/mV) 
  bh = 19.998140294657137     (1) 
  vhh = -87.43164521699035     (mV) 
  Ah = 207.32293925621843     (/ms) 
  b1h = -0.10156699273120488     (/mV) 
  c1h = 0.0009675955808849533     (/mV2) 
  d1h = -2.805630549012438e-06     (/mV3) 
  b2h = -0.15406132188567467     (/mV) 
  c2h = -0.014602550986090954     (/mV2) 
  d2h = -0.0005836270449390211     (/mV3) 

  am = 0.16127412962727178     (/mV) 
  bm = -9.031591516551037     (1) 
  vhm = -65.96424341000792     (mV) 
  Am = 8.387313448352812     (/ms) 
  b1m = 0.02560144405664537     (/mV) 
  c1m = -0.0003394950621593115     (/mV2) 
  d1m = 1.594797554675474e-06     (/mV3) 
  b2m = 0.07623126212839844     (/mV) 
  c2m = -0.0005132814594217189     (/mV2) 
  d2m = 1.1231911005578292e-06     (/mV3) 
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