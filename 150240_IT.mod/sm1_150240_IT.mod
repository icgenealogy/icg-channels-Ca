NEURON
{
  SUFFIX it 
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

  ah = -0.24996366155582864     (/mV) 
  bh = 20.74715829976017     (1) 
  vhh = -91.03113039504531     (mV) 
  Ah = 122.91457215180556     (/ms) 
  b1h = -0.09572599274191393     (/mV) 
  c1h = 0.0009194433204923795     (/mV2) 
  d1h = -2.6566661543123385e-06     (/mV3) 
  b2h = -0.4535503480371976     (/mV) 
  c2h = -0.15519914124640583     (/mV2) 
  d2h = -0.012297808831327427     (/mV3) 
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
}

STATE
{
  h
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*h
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  h' = (hInf - h) / hTau 
}

INITIAL
{
  rates(v)
  h = hInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 


  UNITSON
}