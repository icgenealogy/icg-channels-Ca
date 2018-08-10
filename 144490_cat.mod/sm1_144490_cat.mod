NEURON
{
  SUFFIX cat 
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

  ah = -0.1518481777046413     (/mV) 
  bh = 10.194526812279276     (1) 
  vhh = -68.69579442621071     (mV) 
  Ah = 5571.68598900847     (/ms) 
  b1h = -0.08962072593406151     (/mV) 
  c1h = -0.00030513168810269786     (/mV2) 
  d1h = 3.3458130045808198e-06     (/mV3) 
  b2h = -0.06081674574010138     (/mV) 
  c2h = -0.00022825451910449523     (/mV2) 
  d2h = -2.4016930225400687e-06     (/mV3) 

  am = 0.12106763032539963     (/mV) 
  bm = -2.692152255118588     (1) 
  vhm = -18.923063417220526     (mV) 
  Am = 7.520940734523824     (/ms) 
  b1m = 0.031030968356359302     (/mV) 
  c1m = -0.00024476207323638553     (/mV2) 
  d1m = -1.549853795304591e-06     (/mV3) 
  b2m = 0.08846146713749657     (/mV) 
  c2m = -0.0005653981609059239     (/mV2) 
  d2m = 1.2145304925881319e-06     (/mV3) 
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
  g = gbar*h*h*m*m
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