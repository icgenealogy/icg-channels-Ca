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

  ah = -0.15184748049330596     (/mV) 
  bh = 10.195704383676768     (1) 
  vhh = -68.50664162778986     (mV) 
  Ah = 5549.618134240463     (/ms) 
  b1h = 0.059560842100968535     (/mV) 
  c1h = 0.00018163326894228248     (/mV2) 
  d1h = 1.7813087267142183e-06     (/mV3) 
  b2h = 0.09061878451149129     (/mV) 
  c2h = 0.0002877018053233996     (/mV2) 
  d2h = -3.2814788776440675e-06     (/mV3) 

  am = 0.12105735139964896     (/mV) 
  bm = -2.6925284447367885     (1) 
  vhm = -18.893437569419824     (mV) 
  Am = 7.50196653321117     (/ms) 
  b1m = 0.030447802082116334     (/mV) 
  c1m = -0.0002645316296100796     (/mV2) 
  d1m = -1.7229861113263976e-06     (/mV3) 
  b2m = 0.08793783297533507     (/mV) 
  c2m = -0.0005476277503626256     (/mV2) 
  d2m = 1.0591849889699704e-06     (/mV3) 
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