NEURON
{
  SUFFIX cacurrent 
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

  as = 0.12862304221848012     (/mV) 
  bs = -2.553820985287716     (1) 
  vhs = -15.958559662812448     (mV) 
  As = 4.176449680310365     (/ms) 
  b1s = 0.056602731379993144     (/mV) 
  c1s = 0.0006917170541907877     (/mV2) 
  d1s = 3.6830133732446657e-06     (/mV3) 
  b2s = 0.07111589924530023     (/mV) 
  c2s = -0.0008789707282497681     (/mV2) 
  d2s = 3.536974745624303e-06     (/mV3) 
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
}

STATE
{
  s
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*s
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  s' = (sInf - s) / sTau 
}

INITIAL
{
  rates(v)
  s = sInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 


  UNITSON
}