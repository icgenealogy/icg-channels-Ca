NEURON
{
  SUFFIX car 
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

  ah = -0.0689654553773179     (/mV) 
  bh = 5.427584129954027     (1) 
  vhh = -34.823008446684746     (mV) 
  Ah = 23.336519201800257     (/ms) 
  b1h = -0.04048117069539993     (/mV) 
  c1h = 0.00024893392406527055     (/mV2) 
  d1h = -2.707148124607146e-07     (/mV3) 
  b2h = -0.06694885431457356     (/mV) 
  c2h = -0.00026837755600379567     (/mV2) 
  d2h = 5.7333794597441855e-06     (/mV3) 

  am = 0.17241339562144042     (/mV) 
  bm = -2.586188835044401     (1) 
  vhm = -40.0     (mV) 
  Am = 23.979999956022105     (/ms) 
  b1m = -19.94547018676614     (/mV) 
  c1m = -0.0023177061590933194     (/mV2) 
  d1m = 1.2996760434680711e-05     (/mV3) 
  b2m = -19.877005776614432     (/mV) 
  c2m = 0.0011709702203552164     (/mV2) 
  d2m = -1.3213068228205453e-05     (/mV3) 
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
  g = gbar*h*m
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