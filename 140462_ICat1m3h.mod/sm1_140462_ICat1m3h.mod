NEURON
{
  SUFFIX iCat1m3h 
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

  ah = -0.16554268145539686     (/mV) 
  bh = 12.155975027338775     (1) 
  vhh = -20.39426894954843     (mV) 
  Ah = 66.56209547754645     (/ms) 
  b1h = -0.008491668016427206     (/mV) 
  c1h = -6.161527295623838e-05     (/mV2) 
  d1h = 8.890184498469052e-07     (/mV3) 
  b2h = 0.00848809683949897     (/mV) 
  c2h = -0.0017877781543263755     (/mV2) 
  d2h = 1.3767826638026219e-05     (/mV3) 

  am = 0.1383124035655886     (/mV) 
  bm = -7.114788434277689     (1) 
  vhm = -9.793321035078263     (mV) 
  Am = 4.452091706232902     (/ms) 
  b1m = -0.007781061006168403     (/mV) 
  c1m = 0.0005708542186330583     (/mV2) 
  d1m = 4.536645839792039e-06     (/mV3) 
  b2m = 0.00778303522879308     (/mV) 
  c2m = 8.349959054802286e-05     (/mV2) 
  d2m = -9.302074215035749e-07     (/mV3) 
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
  g = gbar*h*m*m*m
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