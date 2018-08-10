NEURON
{
  SUFFIX gfbp 
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

  ac = 0.1700729369503765     (/mV) 
  bc = -3.6391820189393527     (1) 
  vhc = -24.06004484994703     (mV) 
  Ac = 0.5829159125215874     (/ms) 
  b1c = 0.14025194577105127     (/mV) 
  c1c = 0.0017920058456689714     (/mV2) 
  d1c = 3.103371801489988e-05     (/mV3) 
  b2c = 0.04873997119394232     (/mV) 
  c2c = -0.00028855028373739593     (/mV2) 
  d2c = 5.877615578460363e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  cInf 
  cTau 
}

STATE
{
  c
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*c*c*c
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  c' = (cInf - c) / cTau 
}

INITIAL
{
  rates(v)
  c = cInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    cInf = 1/(1 + exp(-ac*v + bc)) 
    cTau = Ac / ( exp(-(b1c*(v-vhc) + c1c*(v-vhc)^2 + d1c*(v-vhc)^3)) + exp((b2c*(v-vhc) + c2c*(v-vhc)^2 + d2c*(v-vhc)^3)) ) 


  UNITSON
}