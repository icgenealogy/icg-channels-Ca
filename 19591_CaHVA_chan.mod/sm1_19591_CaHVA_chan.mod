NEURON
{
  SUFFIX cahva_chan 
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

  ah = -0.10837949922550341     (/mV) 
  bh = 4.681360746880843     (1) 
  vhh = -71.25928937682262     (mV) 
  Ah = 400.0213681549011     (/ms) 
  b1h = -0.0026097363880455385     (/mV) 
  c1h = 5.226883762395704e-06     (/mV2) 
  d1h = -4.669416150506174e-09     (/mV3) 
  b2h = -0.002609934888393942     (/mV) 
  c2h = -1.5645596874702898e-06     (/mV2) 
  d2h = 5.184217164842922e-09     (/mV3) 

  am = 0.1286160588258156     (/mV) 
  bm = -2.553597960556579     (1) 
  vhm = -15.958559662812448     (mV) 
  Am = 4.176449680310365     (/ms) 
  b1m = 0.056602731379993144     (/mV) 
  c1m = 0.0006917170541907877     (/mV2) 
  d1m = 3.6830133732446657e-06     (/mV3) 
  b2m = 0.07111589924530023     (/mV) 
  c2m = -0.0008789707282497681     (/mV2) 
  d2m = 3.536974745624303e-06     (/mV3) 
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