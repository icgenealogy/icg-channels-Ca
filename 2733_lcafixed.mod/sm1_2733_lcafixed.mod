NEURON
{
  SUFFIX lcafixed 
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

  as = 0.1447679709004134     (/mV) 
  bs = 0.4554304241016674     (1) 
  vhs = -173.94876453379106     (mV) 
  As = 0.023874920939284887     (/ms) 
  b1s = 0.11990846502596438     (/mV) 
  c1s = -0.0013452905534690265     (/mV2) 
  d1s = 4.9775648258015585e-06     (/mV3) 
  b2s = -0.11992045583777129     (/mV) 
  c2s = 0.0009202052490347864     (/mV2) 
  d2s = -1.8597483610240417e-06     (/mV3) 

  ar = -0.12358017152413583     (/mV) 
  br = 3.832466063933877     (1) 
  vhr = -23.568927853822032     (mV) 
  Ar = 225.17107748324335     (/ms) 
  b1r = 0.020624314391063256     (/mV) 
  c1r = 0.0003434688146719652     (/mV2) 
  d1r = 1.9224672406360285e-06     (/mV3) 
  b2r = 0.10623896712773441     (/mV) 
  c2r = -0.0013666124829993164     (/mV2) 
  d2r = 5.561182083934832e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  sInf 
  sTau 
  rInf 
  rTau 
}

STATE
{
  s
  r
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*s*r
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  s' = (sInf - s) / sTau 
  r' = (rInf - r) / rTau 
}

INITIAL
{
  rates(v)
  s = sInf 
  r = rInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 

    rInf = 1/(1 + exp(-ar*v + br)) 
    rTau = Ar / ( exp(-(b1r*(v-vhr) + c1r*(v-vhr)^2 + d1r*(v-vhr)^3)) + exp((b2r*(v-vhr) + c2r*(v-vhr)^2 + d2r*(v-vhr)^3)) ) 


  UNITSON
}