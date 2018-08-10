NEURON
{
  SUFFIX HVA 
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

  au = 0.08849550180154095     (/mV) 
  bu = -2.1769849261611594     (1) 
  vhu = -49.8580815065054     (mV) 
  Au = 2.3356697040315693     (/ms) 
  b1u = -0.01685935720850734     (/mV) 
  c1u = -0.0001493369742641708     (/mV2) 
  d1u = 5.956880015115196e-07     (/mV3) 
  b2u = -0.040107713330717255     (/mV) 
  c2u = -0.00014028474122618006     (/mV2) 
  d2u = -8.519125170173945e-07     (/mV3) 

  az = -0.052890123748525655     (/mV) 
  bz = 0.6666903183620361     (1) 
  vhz = -41.78898823498683     (mV) 
  Az = 280.0200001242232     (/ms) 
  b1z = 1.2276791097971697e-09     (/mV) 
  c1z = -6.536617755574119e-09     (/mV2) 
  d1z = -3.568161593354749e-11     (/mV3) 
  b2z = 1.300062007336963e-09     (/mV) 
  c2z = -6.535812707451669e-09     (/mV2) 
  d2z = -3.570806000596278e-11     (/mV3) 
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