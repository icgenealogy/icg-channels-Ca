NEURON
{
  SUFFIX nca 
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

  ac = 0.12491989001502267     (/mV) 
  bc = -2.6508789426250927     (1) 
  vhc = 19.519151186831856     (mV) 
  Ac = 0.128045885269276     (/ms) 
  b1c = -0.014175743373298721     (/mV) 
  c1c = -0.00020443820056929407     (/mV2) 
  d1c = 1.6169818091294113e-06     (/mV3) 
  b2c = -0.2703193218129928     (/mV) 
  c2c = -0.013485104384535606     (/mV2) 
  d2c = -0.00010219012643843083     (/mV3) 

  ad = -0.12060002118971683     (/mV) 
  bd = 4.83759835973477     (1) 
  vhd = -42.94691416163987     (mV) 
  Ad = 102.91609982056534     (/ms) 
  b1d = 0.020068485753431856     (/mV) 
  c1d = -0.00018418998593998685     (/mV2) 
  d1d = -3.2819004866765165e-06     (/mV3) 
  b2d = 0.07852912745329917     (/mV) 
  c2d = 0.0007753778909264052     (/mV2) 
  d2d = -8.540884494534021e-06     (/mV3) 
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