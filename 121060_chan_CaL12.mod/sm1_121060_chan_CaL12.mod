NEURON
{
  SUFFIX CaL12 
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

  ah = -0.08403177886843627     (/mV) 
  bh = 1.1260370079263922     (1) 
  vhh = -60.77617980703009     (mV) 
  Ah = 29.56000001903031     (/ms) 
  b1h = -9.003551546476147e-08     (/mV) 
  c1h = 1.478673798636503e-09     (/mV2) 
  d1h = -5.119272435946441e-11     (/mV3) 
  b2h = -9.003027784183099e-08     (/mV) 
  c2h = 1.4803463707821689e-09     (/mV2) 
  d2h = -5.121071666567746e-11     (/mV3) 

  am = 0.14925347099477665     (/mV) 
  bm = -1.3283430089439865     (1) 
  vhm = 3.665308986583909     (mV) 
  Am = 0.5306445153853682     (/ms) 
  b1m = -0.058410889303705645     (/mV) 
  c1m = 0.0004640532770890054     (/mV2) 
  d1m = -2.3581308372780397e-06     (/mV3) 
  b2m = -0.02208870455985766     (/mV) 
  c2m = -2.4724271022453868e-05     (/mV2) 
  d2m = 3.0722322260087925e-07     (/mV3) 
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