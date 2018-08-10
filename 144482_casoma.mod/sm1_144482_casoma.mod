NEURON
{
  SUFFIX casoma 
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

  ah = -0.11570324551534167     (/mV) 
  bh = 3.770816581327207     (1) 
  vhh = -33.79235022768801     (mV) 
  Ah = 815.5316358938785     (/ms) 
  b1h = -0.0999475228197242     (/mV) 
  c1h = 0.0013090123749622058     (/mV2) 
  d1h = -5.224427196742729e-06     (/mV3) 
  b2h = -0.02024382998382428     (/mV) 
  c2h = 3.122028309389093e-05     (/mV2) 
  d2h = 5.474466477830297e-07     (/mV3) 

  am = 0.24929600088977238     (/mV) 
  bm = -13.6548256917174     (1) 
  vhm = -55.0     (mV) 
  Am = 37.32000022036888     (/ms) 
  b1m = -12.063846676150423     (/mV) 
  c1m = -0.002630433649100547     (/mV2) 
  d1m = 1.0871652571647767e-05     (/mV3) 
  b2m = -12.083930061207155     (/mV) 
  c2m = 0.0011415488613175657     (/mV2) 
  d2m = -1.1015374132869449e-05     (/mV3) 
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