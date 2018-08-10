NEURON
{
  SUFFIX it2 
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

  ah = -0.19991753387585023     (/mV) 
  bh = 15.594100815061749     (1) 
  vhh = -84.69275005381371     (mV) 
  Ah = 411.163259984288     (/ms) 
  b1h = 0.14579469019892854     (/mV) 
  c1h = 0.011268773870029322     (/mV2) 
  d1h = 0.00035996735689716226     (/mV3) 
  b2h = 0.10722798709276464     (/mV) 
  c2h = -0.001110841408206825     (/mV2) 
  d2h = 3.4007518141588903e-06     (/mV3) 

  am = 0.13513498858959788     (/mV) 
  bm = -6.756746193023375     (1) 
  vhm = -39.248886937524055     (mV) 
  Am = 4.488114815646338     (/ms) 
  b1m = -0.0020966508516964853     (/mV) 
  c1m = -0.00036131049112898627     (/mV2) 
  d1m = 2.313700947614589e-06     (/mV3) 
  b2m = 0.09378182333295686     (/mV) 
  c2m = -0.0019259127562029545     (/mV2) 
  d2m = 9.294869085188484e-06     (/mV3) 
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