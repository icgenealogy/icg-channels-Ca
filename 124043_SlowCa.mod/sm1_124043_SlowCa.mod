NEURON
{
  SUFFIX sca 
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

  ah = -0.04674980090984912     (/mV) 
  bh = 2.690303292713082     (1) 
  vhh = -7.520854447910542     (mV) 
  Ah = 459.5263972886578     (/ms) 
  b1h = -0.012631880865964468     (/mV) 
  c1h = -7.406383821851147e-05     (/mV2) 
  d1h = 8.644459561210136e-07     (/mV3) 
  b2h = 0.012632064191421375     (/mV) 
  c2h = -0.00020226208512719038     (/mV2) 
  d2h = -2.1284215647998834e-06     (/mV3) 

  am = 0.21843414186295607     (/mV) 
  bm = -7.216315114895137     (1) 
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