NEURON
{
  SUFFIX cahi 
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

  ax = 0.16000023019079943     (/mV) 
  bx = -4.800002420619134     (1) 
  vhx = -73.56183466673515     (mV) 
  Ax = 10.020000007078417     (/ms) 
  b1x = -3.322440488619548e-09     (/mV) 
  c1x = 1.803021258416065e-09     (/mV2) 
  d1x = -5.310905572541754e-11     (/mV3) 
  b2x = -3.3883273721517777e-09     (/mV) 
  c2x = 1.8063835814949279e-09     (/mV2) 
  d2x = -5.3134583518005927e-11     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  xInf 
  xTau 
}

STATE
{
  x
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*x*x*x
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  x' = (xInf - x) / xTau 
}

INITIAL
{
  rates(v)
  x = xInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    xInf = 1/(1 + exp(-ax*v + bx)) 
    xTau = Ax / ( exp(-(b1x*(v-vhx) + c1x*(v-vhx)^2 + d1x*(v-vhx)^3)) + exp((b2x*(v-vhx) + c2x*(v-vhx)^2 + d2x*(v-vhx)^3)) ) 


  UNITSON
}