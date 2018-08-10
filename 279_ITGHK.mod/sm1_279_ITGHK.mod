NEURON
{
  SUFFIX itGHK 
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

  ah = -0.24996475000279061     (/mV) 
  bh = 20.747259887229156     (1) 
  vhh = -90.67398048797398     (mV) 
  Ah = 138.4627381730996     (/ms) 
  b1h = 0.16491461149836784     (/mV) 
  c1h = 0.024890350566865463     (/mV2) 
  d1h = 0.00155708032946945     (/mV3) 
  b2h = 0.09597638172760893     (/mV) 
  c2h = -0.0008901294539839759     (/mV2) 
  d2h = 2.5189302379913047e-06     (/mV3) 

  am = 0.16129001805188078     (/mV) 
  bm = -9.516109813410614     (1) 
  vhm = -66.88781736827151     (mV) 
  Am = 2.943039143293572     (/ms) 
  b1m = -0.07734211185294798     (/mV) 
  c1m = 0.0005501849487859239     (/mV2) 
  d1m = -1.295000174295233e-06     (/mV3) 
  b2m = -0.01785459484051065     (/mV) 
  c2m = 0.0003999317574392579     (/mV2) 
  d2m = -2.642794699461656e-06     (/mV3) 
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