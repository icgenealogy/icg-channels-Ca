NEURON
{
  SUFFIX can 
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

  ah = -0.1144890433473227     (/mV) 
  bh = 7.7346722839377415     (1) 
  vhh = -54.09577936784464     (mV) 
  Ah = 27.169016402446804     (/ms) 
  b1h = 0.017361279819092558     (/mV) 
  c1h = -0.0005202670696740075     (/mV2) 
  d1h = 2.4779910994697635e-06     (/mV3) 
  b2h = 0.05838613013550604     (/mV) 
  c2h = -0.0007824499833314274     (/mV2) 
  d2h = 2.9751596197267463e-06     (/mV3) 

  am = 0.16622380812229715     (/mV) 
  bm = -0.6110630412556485     (1) 
  vhm = -9.268933283973142     (mV) 
  Am = 0.6849545329189091     (/ms) 
  b1m = -0.03714523756798898     (/mV) 
  c1m = 0.000258840134199443     (/mV2) 
  d1m = -8.029193849519826e-07     (/mV3) 
  b2m = -0.1351127774756355     (/mV) 
  c2m = -0.001572885959651985     (/mV2) 
  d2m = -1.8136304933855764e-05     (/mV3) 
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