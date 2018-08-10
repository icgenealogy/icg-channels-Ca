NEURON
{
  SUFFIX ch_CavN 
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

  ac = 0.12491996446477924     (/mV) 
  bc = -2.6508798355136363     (1) 
  vhc = -18.106967913497417     (mV) 
  Ac = 5.23284122106072     (/ms) 
  b1c = 0.03481937944608435     (/mV) 
  c1c = -0.0002717165426902642     (/mV2) 
  d1c = -1.8199704475281875e-06     (/mV3) 
  b2c = 0.08756578042030218     (/mV) 
  c2c = -0.0005570092254326926     (/mV2) 
  d2c = 1.1386705747008837e-06     (/mV3) 

  ad = -0.12027422725228726     (/mV) 
  bd = 4.8274898621560105     (1) 
  vhd = -43.021309590212205     (mV) 
  Ad = 2168.549531667855     (/ms) 
  b1d = 0.027780304528536595     (/mV) 
  c1d = 0.00012358392328154583     (/mV2) 
  d1d = 7.629870239641727e-07     (/mV3) 
  b2d = 0.08612892307855889     (/mV) 
  c2d = 0.00039340504177351145     (/mV2) 
  d2d = -4.300025551587843e-06     (/mV3) 
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
  dInf 
  dTau 
}

STATE
{
  c
  d
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*c*c*d
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  c' = (cInf - c) / cTau 
  d' = (dInf - d) / dTau 
}

INITIAL
{
  rates(v)
  c = cInf 
  d = dInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    cInf = 1/(1 + exp(-ac*v + bc)) 
    cTau = Ac / ( exp(-(b1c*(v-vhc) + c1c*(v-vhc)^2 + d1c*(v-vhc)^3)) + exp((b2c*(v-vhc) + c2c*(v-vhc)^2 + d2c*(v-vhc)^3)) ) 

    dInf = 1/(1 + exp(-ad*v + bd)) 
    dTau = Ad / ( exp(-(b1d*(v-vhd) + c1d*(v-vhd)^2 + d1d*(v-vhd)^3)) + exp((b2d*(v-vhd) + c2d*(v-vhd)^2 + d2d*(v-vhd)^3)) ) 


  UNITSON
}