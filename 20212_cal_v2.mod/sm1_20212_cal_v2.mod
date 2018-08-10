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

  am = 0.16666712993226562     (/mV) 
  bm = -1.49999962931087     (1) 
  vhm = -10.346956048031526     (mV) 
  Am = 0.6723398197035307     (/ms) 
  b1m = -0.03459133916053985     (/mV) 
  c1m = -0.00020730860991972726     (/mV2) 
  d1m = 4.36307629218726e-06     (/mV3) 
  b2m = -0.03838746593703344     (/mV) 
  c2m = 7.267885618296045e-05     (/mV2) 
  d2m = 3.962714196390634e-06     (/mV3) 
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
  g = gbar*m*m
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