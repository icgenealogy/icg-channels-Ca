NEURON
{
  SUFFIX Calcium 
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

  an = 0.08849970845930746     (/mV) 
  bn = -1.5311881646949834     (1) 
  vhn = -65.16259053217236     (mV) 
  An = 6.020000000224573     (/ms) 
  b1n = 1.5563289559202296e-08     (/mV) 
  c1n = -8.705175265401171e-09     (/mV2) 
  d1n = 6.141864735877141e-11     (/mV3) 
  b2n = 1.556002765776949e-08     (/mV) 
  c2n = -8.705145723414579e-09     (/mV2) 
  d2n = 6.14184359838454e-11     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  nInf 
  nTau 
}

STATE
{
  n
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*n
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  n' = (nInf - n) / nTau 
}

INITIAL
{
  rates(v)
  n = nInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    nInf = 1/(1 + exp(-an*v + bn)) 
    nTau = An / ( exp(-(b1n*(v-vhn) + c1n*(v-vhn)^2 + d1n*(v-vhn)^3)) + exp((b2n*(v-vhn) + c2n*(v-vhn)^2 + d2n*(v-vhn)^3)) ) 


  UNITSON
}