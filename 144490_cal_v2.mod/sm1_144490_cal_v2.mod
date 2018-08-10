NEURON
{
  SUFFIX cal 
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

  am = 0.20100939900997114     (/mV) 
  bm = -5.972199966105519     (1) 
  vhm = -35.63632144938852     (mV) 
  Am = 1.6577770607258873     (/ms) 
  b1m = 0.1076470460921503     (/mV) 
  c1m = 0.001677558019947095     (/mV2) 
  d1m = 1.6670453632213226e-05     (/mV3) 
  b2m = 0.1081831216175252     (/mV) 
  c2m = -0.0012084790763570435     (/mV2) 
  d2m = 4.755820069190982e-06     (/mV3) 
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
  g = gbar*m
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