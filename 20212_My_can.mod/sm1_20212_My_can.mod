NEURON
{
  SUFFIX mycan 
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

  ah = -0.10870515972930393     (/mV) 
  bh = 4.783109288517836     (1) 
  vhh = 4.324959560695309     (mV) 
  Ah = 142.62738419944276     (/ms) 
  b1h = -0.05716093273841285     (/mV) 
  c1h = 0.0004967185600703535     (/mV2) 
  d1h = -1.870246002677475e-06     (/mV3) 
  b2h = -0.04615075991434296     (/mV) 
  c2h = 8.503094046559675e-05     (/mV2) 
  d2h = 5.263803004131603e-07     (/mV3) 

  am = 0.12047590513667544     (/mV) 
  bm = -1.3257947220043511     (1) 
  vhm = -21.983750687086243     (mV) 
  Am = 36.527498817543844     (/ms) 
  b1m = -0.06690461958786642     (/mV) 
  c1m = 0.0011906758823871837     (/mV2) 
  d1m = -4.998449538745474e-06     (/mV3) 
  b2m = 0.016477822773787735     (/mV) 
  c2m = -0.0005153384481894171     (/mV2) 
  d2m = -2.7316302848394033e-06     (/mV3) 
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