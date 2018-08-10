NEURON
{
  SUFFIX it2 
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

  ah = -0.15619194182544766     (/mV) 
  bh = 10.88691627701481     (1) 
  vhh = -69.22940834551174     (mV) 
  Ah = 175.84109335014202     (/ms) 
  b1h = -0.12179393444304837     (/mV) 
  c1h = 0.0017443315868462702     (/mV2) 
  d1h = -6.527960603634399e-06     (/mV3) 
  b2h = 0.009274079882419309     (/mV) 
  c2h = 0.00033095926196336693     (/mV2) 
  d2h = -1.905504890914134e-06     (/mV3) 

  am = 0.11494242170455704     (/mV) 
  bm = -4.390795093618215     (1) 
  vhm = -43.57853715739607     (mV) 
  Am = 13.489211994674092     (/ms) 
  b1m = 0.0661724209131739     (/mV) 
  c1m = 0.0008102448228383745     (/mV2) 
  d1m = 4.948997654190019e-06     (/mV3) 
  b2m = 0.07728064531624394     (/mV) 
  c2m = -0.0008000062927314536     (/mV2) 
  d2m = 2.658632734340828e-06     (/mV3) 
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