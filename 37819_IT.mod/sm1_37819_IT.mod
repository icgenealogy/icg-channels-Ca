NEURON
{
  SUFFIX ittc 
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

  ah = -0.2499636615557238     (/mV) 
  bh = 20.74715829975128     (1) 
  vhh = -90.91958565008355     (mV) 
  Ah = 121.58741951990578     (/ms) 
  b1h = -0.09549890061477229     (/mV) 
  c1h = 0.0009184086518089025     (/mV2) 
  d1h = -2.656246131315173e-06     (/mV3) 
  b2h = -0.5356627355730413     (/mV) 
  c2h = -0.1823222014040699     (/mV2) 
  d2h = -0.014176136461621866     (/mV3) 

  am = 0.161290018051897     (/mV) 
  bm = -9.51610981341163     (1) 
  vhm = -67.36099350444366     (mV) 
  Am = 6.0678455170859396     (/ms) 
  b1m = 0.026868501365355648     (/mV) 
  c1m = -0.00037700726114250584     (/mV2) 
  d1m = 6.840824393552443e-07     (/mV3) 
  b2m = 0.07465311522549167     (/mV) 
  c2m = -0.00047323038126702646     (/mV2) 
  d2m = 8.48924674426705e-07     (/mV3) 
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