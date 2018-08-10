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

  ah = -0.2499158872273023     (/mV) 
  bh = 19.993540359649963     (1) 
  vhh = -86.89103721942759     (mV) 
  Ah = 198.4967887782065     (/ms) 
  b1h = -0.10116410892534494     (/mV) 
  c1h = 0.0009706613432781962     (/mV2) 
  d1h = -2.8310104542557655e-06     (/mV3) 
  b2h = -0.15966780314763923     (/mV) 
  c2h = -0.015539120992902972     (/mV2) 
  d2h = -0.0006130777042196039     (/mV3) 

  am = 0.16127681368474653     (/mV) 
  bm = -9.031522495670728     (1) 
  vhm = -65.52472284013075     (mV) 
  Am = 8.292826823865791     (/ms) 
  b1m = -0.0768463014070534     (/mV) 
  c1m = 0.0005198759128190201     (/mV2) 
  d1m = -1.1243100323050152e-06     (/mV3) 
  b2m = -0.024132757837549347     (/mV) 
  c2m = 0.0003502317385487691     (/mV2) 
  d2m = -1.8077225319535331e-06     (/mV3) 
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