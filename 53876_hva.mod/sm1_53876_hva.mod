NEURON
{
  SUFFIX hva 
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

  as = 0.12385149107567446     (/mV) 
  bs = -2.5117584350939297     (1) 
  vhs = -19.08769834046402     (mV) 
  As = 1.8209019715813302     (/ms) 
  b1s = -0.040825251520093156     (/mV) 
  c1s = 0.0004433718152028443     (/mV2) 
  d1s = -1.6088904198758467e-06     (/mV3) 
  b2s = -0.08374376219164005     (/mV) 
  c2s = -0.0012254889586799345     (/mV2) 
  d2s = -7.197037751651765e-06     (/mV3) 

  ar = -0.10840068980482528     (/mV) 
  br = 4.682023298949572     (1) 
  vhr = -356.33637684716376     (mV) 
  Ar = 226.66963122296715     (/ms) 
  b1r = -0.0038872276269353817     (/mV) 
  c1r = 6.722106418851732e-06     (/mV2) 
  d1r = -4.092120754395216e-09     (/mV3) 
  b2r = -0.000166718234455032     (/mV) 
  c2r = -1.7801297063247158e-05     (/mV2) 
  d2r = 1.7688715648627237e-08     (/mV3) 
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
  g = gbar*s*s*r
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