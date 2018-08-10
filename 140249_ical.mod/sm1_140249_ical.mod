NEURON
{
  SUFFIX ical 
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

  am = 0.112219636326872     (/mV) 
  bm = -0.843058718823063     (1) 
  vhm = -19.63618358718949     (mV) 
  Am = 0.5016287851665256     (/ms) 
  b1m = 0.03676006070964409     (/mV) 
  c1m = 4.928892363011117e-05     (/mV2) 
  d1m = -3.6140983965169776e-06     (/mV3) 
  b2m = 0.01832800914030393     (/mV) 
  c2m = 0.00022875323754744701     (/mV2) 
  d2m = -2.773266445497444e-06     (/mV3) 
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