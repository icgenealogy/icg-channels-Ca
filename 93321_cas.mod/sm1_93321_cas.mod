NEURON
{
  SUFFIX cas 
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

  ah = -0.16124442882055207     (/mV) 
  bh = 9.675103560273477     (1) 
  vhh = -40.8340040452369     (mV) 
  Ah = 194.7871019671599     (/ms) 
  b1h = -0.002711739323915341     (/mV) 
  c1h = -0.0002679813866510584     (/mV2) 
  d1h = 1.803102088467332e-06     (/mV3) 
  b2h = 0.06862522345712235     (/mV) 
  c2h = -0.0016772780371139455     (/mV2) 
  d2h = 8.593619639575351e-06     (/mV3) 

  am = 0.12345713118040727     (/mV) 
  bm = -4.074118404578291     (1) 
  vhm = -43.694250819740745     (mV) 
  Am = 46.76490357238783     (/ms) 
  b1m = 0.06801579486433466     (/mV) 
  c1m = 2.15831496364481e-05     (/mV2) 
  d1m = -3.2688673609582098e-06     (/mV3) 
  b2m = 0.12195024731250376     (/mV) 
  c2m = -0.0013223409280260934     (/mV2) 
  d2m = 4.524589183261423e-06     (/mV3) 
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