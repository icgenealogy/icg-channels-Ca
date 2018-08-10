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

  ac = 0.1249231047709968     (/mV) 
  bc = -2.650654659379263     (1) 
  vhc = 19.392509624054433     (mV) 
  Ac = 0.12739962137671443     (/ms) 
  b1c = 0.2731483492727873     (/mV) 
  c1c = 0.013725538684716553     (/mV2) 
  d1c = 0.00010420868755341969     (/mV3) 
  b2c = 0.014109624626170615     (/mV) 
  c2c = 0.00020256964510678344     (/mV2) 
  d2c = -1.6253702582700024e-06     (/mV3) 

  ad = -0.07914129117391691     (/mV) 
  bd = 4.8300657745040025     (1) 
  vhd = -120.58367625545357     (mV) 
  Ad = 2546.5301313719406     (/ms) 
  b1d = 1.0269129299317343     (/mV) 
  c1d = 0.0033431694101741806     (/mV2) 
  d1d = -7.412776538699895e-06     (/mV3) 
  b2d = 0.02298788302768353     (/mV) 
  c2d = -0.00014619238086087416     (/mV2) 
  d2d = 4.97133315046955e-06     (/mV3) 
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