NEURON
{
  SUFFIX cadyn 
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

  au = 0.08849542573689526     (/mV) 
  bu = -2.176982363317632     (1) 
  vhu = -43.64309399763411     (mV) 
  Au = 2.4584208053676817     (/ms) 
  b1u = -0.024629889805879986     (/mV) 
  c1u = -6.401271816101378e-05     (/mV2) 
  d1u = 2.5392147486314133e-07     (/mV3) 
  b2u = -0.03724050355132817     (/mV) 
  c2u = -9.690962132770928e-05     (/mV2) 
  d2u = -6.057522594858922e-07     (/mV3) 

  az = -0.052847444878612926     (/mV) 
  bz = 0.6667499957599244     (1) 
  vhz = 0.39773325585906133     (mV) 
  Az = 840.0040191507953     (/ms) 
  b1z = -0.003349878230929574     (/mV) 
  c1z = 6.001474936918206e-06     (/mV2) 
  d1z = -4.444889274384299e-09     (/mV3) 
  b2z = -0.0033491070740590367     (/mV) 
  c2z = -5.22410393576491e-06     (/mV2) 
  d2z = -1.963103835291043e-09     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  uInf 
  uTau 
  zInf 
  zTau 
}

STATE
{
  u
  z
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*u*u*z
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  u' = (uInf - u) / uTau 
  z' = (zInf - z) / zTau 
}

INITIAL
{
  rates(v)
  u = uInf 
  z = zInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    uInf = 1/(1 + exp(-au*v + bu)) 
    uTau = Au / ( exp(-(b1u*(v-vhu) + c1u*(v-vhu)^2 + d1u*(v-vhu)^3)) + exp((b2u*(v-vhu) + c2u*(v-vhu)^2 + d2u*(v-vhu)^3)) ) 

    zInf = 1/(1 + exp(-az*v + bz)) 
    zTau = Az / ( exp(-(b1z*(v-vhz) + c1z*(v-vhz)^2 + d1z*(v-vhz)^3)) + exp((b2z*(v-vhz) + c2z*(v-vhz)^2 + d2z*(v-vhz)^3)) ) 


  UNITSON
}