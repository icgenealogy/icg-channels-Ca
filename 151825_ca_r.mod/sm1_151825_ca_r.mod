NEURON
{
  SUFFIX car 
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

  ah = -0.1280797679868071     (/mV) 
  bh = 10.1201988213117     (1) 
  vhh = -79.0673902413373     (mV) 
  Ah = 1541.9992188782046     (/ms) 
  b1h = -0.046882809869755335     (/mV) 
  c1h = -1.5140787150753061e-06     (/mV2) 
  d1h = 6.9145611331724854e-09     (/mV3) 
  b2h = -0.04718121233333167     (/mV) 
  c2h = -4.60809356841248e-06     (/mV2) 
  d2h = -5.3991239075362206e-08     (/mV3) 

  am = 0.13513492506969835     (/mV) 
  bm = -3.37836427750202     (1) 
  vhm = -26.198741342628676     (mV) 
  Am = 11.0054408919676     (/ms) 
  b1m = -0.03078103308733034     (/mV) 
  c1m = -1.469058156625617e-05     (/mV2) 
  d1m = 7.267489206532813e-08     (/mV3) 
  b2m = -0.033238680868603813     (/mV) 
  c2m = -1.649463407717889e-05     (/mV2) 
  d2m = -8.314308217191418e-08     (/mV3) 
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