NEURON
{
  SUFFIX ca 
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

  ah = -0.04675941224640286     (/mV) 
  bh = 2.690479599738953     (1) 
  vhh = -7.523230726434089     (mV) 
  Ah = 274.420745263432     (/ms) 
  b1h = -0.012632958830894544     (/mV) 
  c1h = 0.0002022326132042608     (/mV2) 
  d1h = 2.128591396025784e-06     (/mV3) 
  b2h = 0.012631640502296078     (/mV) 
  c2h = 7.406450113661833e-05     (/mV2) 
  d2h = -8.644039198673677e-07     (/mV3) 

  am = 0.21843415346818706     (/mV) 
  bm = -7.21631554301052     (1) 
  vhm = -36.52624248951748     (mV) 
  Am = 7.668928484436402     (/ms) 
  b1m = 0.09313902509434949     (/mV) 
  c1m = 0.0011951189367805707     (/mV2) 
  d1m = 1.2278433874662704e-05     (/mV3) 
  b2m = 0.13811971459975708     (/mV) 
  c2m = -0.0017033553328188909     (/mV2) 
  d2m = 7.196628767444811e-06     (/mV3) 
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