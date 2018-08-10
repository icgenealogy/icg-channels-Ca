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

  ah = -0.24811855758028756     (/mV) 
  bh = 15.8796925730355     (1) 
  vhh = -59.17514464405529     (mV) 
  Ah = 80.00877175815491     (/ms) 
  b1h = 0.009830364041732438     (/mV) 
  c1h = 4.829156531749825e-05     (/mV2) 
  d1h = 8.759511929834821e-08     (/mV3) 
  b2h = 0.009845966400275004     (/mV) 
  c2h = -4.8741697290902874e-05     (/mV2) 
  d2h = 8.403565643269648e-08     (/mV3) 

  am = 0.16129003482236465     (/mV) 
  bm = -6.532240615256907     (1) 
  vhm = -52.52074465887092     (mV) 
  Am = 5.010577495715161     (/ms) 
  b1m = -0.0668455693584409     (/mV) 
  c1m = 0.0003061394059939851     (/mV2) 
  d1m = -7.555766699514173e-08     (/mV3) 
  b2m = -0.04668523065089739     (/mV) 
  c2m = 0.00015262091127039311     (/mV2) 
  d2m = 1.790500421921394e-06     (/mV3) 
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