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

  ah = -0.24996366155646332     (/mV) 
  bh = 20.74715829981285     (1) 
  vhh = -90.91958565008355     (mV) 
  Ah = 121.58741951990578     (/ms) 
  b1h = -0.09549890061477229     (/mV) 
  c1h = 0.0009184086518089025     (/mV2) 
  d1h = -2.656246131315173e-06     (/mV3) 
  b2h = -0.5356627355730413     (/mV) 
  c2h = -0.1823222014040699     (/mV2) 
  d2h = -0.014176136461621866     (/mV3) 
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