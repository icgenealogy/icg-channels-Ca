NEURON
{
  SUFFIX L_HVA_Ca 
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

  am = 0.16666250134549782     (/mV) 
  bm = -1.6666229439384443     (1) 
  vhm = -58.23077533622195     (mV) 
  Am = 40.02000006144667     (/ms) 
  b1m = 4.739094121080306e-08     (/mV) 
  c1m = 6.972353346501731e-09     (/mV2) 
  d1m = 4.1569376829822156e-11     (/mV3) 
  b2m = 4.7440036938700294e-08     (/mV) 
  c2m = 6.976058501424149e-09     (/mV2) 
  d2m = 4.1516397994944993e-11     (/mV3) 
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