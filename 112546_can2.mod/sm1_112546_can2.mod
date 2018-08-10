NEURON
{
  SUFFIX can 
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

  ah = -0.12057609545916702     (/mV) 
  bh = 4.836555896355515     (1) 
  vhh = -23.877883441565054     (mV) 
  Ah = 160.00071379217124     (/ms) 
  b1h = 0.010647642378553809     (/mV) 
  c1h = 5.6183378777828846e-05     (/mV2) 
  d1h = 9.829913265759383e-08     (/mV3) 
  b2h = 0.0106535316656569     (/mV) 
  c2h = -5.7452249248771164e-05     (/mV2) 
  d2h = 1.0887306265427356e-07     (/mV3) 

  am = 0.12503908348420853     (/mV) 
  bm = -2.687513112025324     (1) 
  vhm = -84.87239409524733     (mV) 
  Am = 2.288656970414251     (/ms) 
  b1m = 0.05761567551220665     (/mV) 
  c1m = 0.0003079685360984481     (/mV2) 
  d1m = 5.358154861211142e-06     (/mV3) 
  b2m = 0.015838189837279297     (/mV) 
  c2m = 0.0001564035779015916     (/mV2) 
  d2m = -5.88265751308107e-07     (/mV3) 
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