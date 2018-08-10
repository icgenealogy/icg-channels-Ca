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

  ah = -0.12012487437446057     (/mV) 
  bh = 4.824128848324177     (1) 
  vhh = -43.13906323941121     (mV) 
  Ah = 3023.1809186111086     (/ms) 
  b1h = 0.028034260224344992     (/mV) 
  c1h = 0.00012880658117771807     (/mV2) 
  d1h = 7.998470131002197e-07     (/mV3) 
  b2h = 0.08564576877244418     (/mV) 
  c2h = 0.00040223660932815064     (/mV2) 
  d2h = -4.334210307374408e-06     (/mV3) 

  am = 0.12502800386206978     (/mV) 
  bm = -2.6878072345084574     (1) 
  vhm = -18.224896156445734     (mV) 
  Am = 7.144146302966295     (/ms) 
  b1m = 0.034222898664970004     (/mV) 
  c1m = -0.00027818416179763535     (/mV2) 
  d1m = -1.8421670178989296e-06     (/mV3) 
  b2m = 0.0880166454479708     (/mV) 
  c2m = -0.0005620507614297721     (/mV2) 
  d2m = 1.1597133616210079e-06     (/mV3) 
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