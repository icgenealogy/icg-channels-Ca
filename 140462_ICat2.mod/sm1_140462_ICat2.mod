NEURON
{
  SUFFIX iCat2 
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

  am = 0.11566815147002105     (/mV) 
  bm = -3.9878736492296554     (1) 
  vhm = -33.22722959232927     (mV) 
  Am = 4.538387298004056     (/ms) 
  b1m = 0.06933060565545457     (/mV) 
  c1m = 0.0010884288447066545     (/mV2) 
  d1m = 7.330176474739893e-06     (/mV3) 
  b2m = 0.05061983839134087     (/mV) 
  c2m = -0.0005329531190122065     (/mV2) 
  d2m = 1.8165034056366559e-06     (/mV3) 
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