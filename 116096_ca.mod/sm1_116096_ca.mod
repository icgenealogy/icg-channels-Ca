NEURON
{
  SUFFIX ca 
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

  ah = -0.046765869594829916     (/mV) 
  bh = 2.6906001070746357     (1) 
  vhh = -7.521977928489206     (mV) 
  Ah = 153.1822625147669     (/ms) 
  b1h = -0.01263168436298006     (/mV) 
  c1h = 0.000202249346241891     (/mV2) 
  d1h = 2.1288041047170544e-06     (/mV3) 
  b2h = 0.012631877012663884     (/mV) 
  c2h = 7.405909086724221e-05     (/mV2) 
  d2h = -8.643812862646538e-07     (/mV3) 

  am = 0.21843407363927395     (/mV) 
  bm = -7.21631259999547     (1) 
  vhm = -36.601656586118494     (mV) 
  Am = 4.282852535389412     (/ms) 
  b1m = 0.09315971496619026     (/mV) 
  c1m = 0.0011987421096272713     (/mV2) 
  d1m = 1.2261383194933285e-05     (/mV3) 
  b2m = 0.1371479414654808     (/mV) 
  c2m = -0.0016781137075827695     (/mV2) 
  d2m = 7.017224621733026e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hInf 
  hTau 
  mInf 
  mTau 
}

STATE
{
  h
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h*m*m
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}