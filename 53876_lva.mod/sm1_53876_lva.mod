NEURON
{
  SUFFIX lva 
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

  ah = -0.15873087134384395     (/mV) 
  bh = 13.25416236874282     (1) 
  vhh = -38.05470903593286     (mV) 
  Ah = 378.8044416341424     (/ms) 
  b1h = -0.03284204611047795     (/mV) 
  c1h = 0.0003249881864749162     (/mV2) 
  d1h = -8.431076716420681e-06     (/mV3) 
  b2h = 0.03281870688746199     (/mV) 
  c2h = 0.0004409806706414635     (/mV2) 
  d2h = -8.026649109999003e-07     (/mV3) 

  am = 0.12820497886691284     (/mV) 
  bm = -8.076913512890354     (1) 
  vhm = -48.25259531983375     (mV) 
  Am = 3.1673808226950504     (/ms) 
  b1m = -0.07539711744480421     (/mV) 
  c1m = 0.001176664945777755     (/mV2) 
  d1m = -4.912473329195458e-06     (/mV3) 
  b2m = -0.006736728060346503     (/mV) 
  c2m = 0.0004623506908603561     (/mV2) 
  d2m = -2.525568295776916e-06     (/mV3) 

  ad = 0.13469472907007218     (/mV) 
  bd = -9.886764427701602     (1) 
  vhd = -68.88476125614604     (mV) 
  Ad = 37.87658029063554     (/ms) 
  b1d = 0.015132618292215494     (/mV) 
  c1d = -0.0006337083317046579     (/mV2) 
  d1d = 2.7472940766016094e-06     (/mV3) 
  b2d = 0.03858291271137853     (/mV) 
  c2d = -0.00027240266278128     (/mV2) 
  d2d = 1.1936620542853077e-06     (/mV3) 
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
  dInf 
  dTau 
}

STATE
{
  h
  m
  d
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
  d' = (dInf - d) / dTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
  d = dInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    dInf = 1/(1 + exp(-ad*v + bd)) 
    dTau = Ad / ( exp(-(b1d*(v-vhd) + c1d*(v-vhd)^2 + d1d*(v-vhd)^3)) + exp((b2d*(v-vhd) + c2d*(v-vhd)^2 + d2d*(v-vhd)^3)) ) 


  UNITSON
}