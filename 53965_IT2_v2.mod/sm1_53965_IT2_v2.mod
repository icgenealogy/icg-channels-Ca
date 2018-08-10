NEURON
{
  SUFFIX it2 
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

  ah = -0.11735165238881919     (/mV) 
  bh = 9.769884547018668     (1) 
  vhh = -92.0915749615856     (mV) 
  Ah = 1640.7631813508212     (/ms) 
  b1h = -0.17584135696867575     (/mV) 
  c1h = 0.0017654665729545651     (/mV2) 
  d1h = -5.241769373358544e-06     (/mV3) 
  b2h = -0.14221474664792766     (/mV) 
  c2h = -0.022919450102962063     (/mV2) 
  d2h = -0.0018840737433364885     (/mV3) 

  am = 0.1273883531382613     (/mV) 
  bm = -6.687886507679375     (1) 
  vhm = -75.61399579338772     (mV) 
  Am = 3.1810044308646854     (/ms) 
  b1m = 0.06424913445611911     (/mV) 
  c1m = 0.0014925696998000944     (/mV2) 
  d1m = 2.4133469163197913e-05     (/mV3) 
  b2m = 0.06181417415823518     (/mV) 
  c2m = -0.0005139612069928713     (/mV2) 
  d2m = 1.3773310250124543e-06     (/mV3) 
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