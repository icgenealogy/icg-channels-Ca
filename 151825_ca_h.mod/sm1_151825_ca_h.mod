NEURON
{
  SUFFIX cah 
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

  ah = -0.047597764272821405     (/mV) 
  bh = 1.7579853243093042     (1) 
  vhh = 19.567649590554907     (mV) 
  Ah = 56.865699785606346     (/ms) 
  b1h = 0.04730460797258119     (/mV) 
  c1h = 3.999991329251671e-06     (/mV2) 
  d1h = 8.94893300327616e-09     (/mV3) 
  b2h = 0.046583044871013946     (/mV) 
  c2h = 8.157538963405303e-06     (/mV2) 
  d2h = -5.7535878441656e-08     (/mV3) 

  am = 0.12232953393001068     (/mV) 
  bm = -2.2566019315140076     (1) 
  vhm = -10.630833023487204     (mV) 
  Am = 0.28347985483320054     (/ms) 
  b1m = -0.04208802284942046     (/mV) 
  c1m = 0.000153434826441664     (/mV2) 
  d1m = -3.6840992151055814e-07     (/mV3) 
  b2m = -0.011425057177964439     (/mV) 
  c2m = 0.00019496457350285133     (/mV2) 
  d2m = 7.651806404516868e-07     (/mV3) 
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