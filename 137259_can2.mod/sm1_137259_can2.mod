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

  ah = -0.12057609545895989     (/mV) 
  bh = 4.836555896346173     (1) 
  vhh = -25.160528963885078     (mV) 
  Ah = 160.0011135835833     (/ms) 
  b1h = -0.010978421265833307     (/mV) 
  c1h = 6.0988539735708146e-05     (/mV2) 
  d1h = -1.1898729349119166e-07     (/mV3) 
  b2h = -0.010972629868861643     (/mV) 
  c2h = -5.969101392077955e-05     (/mV2) 
  d2h = -1.0785407663599841e-07     (/mV3) 

  am = 0.12503908619560147     (/mV) 
  bm = -2.6875131933822964     (1) 
  vhm = -13.531748709793382     (mV) 
  Am = 4.775552488851579     (/ms) 
  b1m = -0.06321030639598807     (/mV) 
  c1m = -0.00023830671926875472     (/mV2) 
  d1m = 3.411080193545486e-06     (/mV3) 
  b2m = -0.0019530738439555205     (/mV) 
  c2m = 0.0001300063746082833     (/mV2) 
  d2m = 9.726488464352131e-07     (/mV3) 
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