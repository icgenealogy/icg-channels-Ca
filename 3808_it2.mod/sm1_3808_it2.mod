NEURON
{
  SUFFIX iT2 
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

  ah = -0.19993542210402634     (/mV) 
  bh = 15.995250200074075     (1) 
  vhh = -73.28717253586484     (mV) 
  Ah = 194.93932200117035     (/ms) 
  b1h = -0.15654944621710706     (/mV) 
  c1h = 0.0026979695981573593     (/mV2) 
  d1h = -1.0777840372452132e-05     (/mV3) 
  b2h = -0.0014644850155499012     (/mV) 
  c2h = 0.0004753230109907685     (/mV2) 
  d2h = -2.4410631341834975e-06     (/mV3) 

  am = 0.13513489289273048     (/mV) 
  bm = -7.027010879107175     (1) 
  vhm = -41.95326371687706     (mV) 
  Am = 1.8820381345497732     (/ms) 
  b1m = -0.09349005460395643     (/mV) 
  c1m = 0.0018741439487199793     (/mV2) 
  d1m = -8.879055271032308e-06     (/mV3) 
  b2m = 0.0010385515093895015     (/mV) 
  c2m = 0.00037546480954605595     (/mV2) 
  d2m = -2.3421247172261605e-06     (/mV3) 
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