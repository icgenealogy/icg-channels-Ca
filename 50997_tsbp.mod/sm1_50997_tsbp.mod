NEURON
{
  SUFFIX tsbp 
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

  ac = 0.13843598156786474     (/mV) 
  bc = -5.2434308036915445     (1) 
  vhc = -42.416444359880046     (mV) 
  Ac = 0.09442609842410497     (/ms) 
  b1c = -0.011323891956057932     (/mV) 
  c1c = -7.343490331246781e-05     (/mV2) 
  d1c = 7.523761808034001e-07     (/mV3) 
  b2c = -0.11974868087558861     (/mV) 
  c2c = 0.0004121901716282639     (/mV2) 
  d2c = 3.228328652680067e-06     (/mV3) 
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