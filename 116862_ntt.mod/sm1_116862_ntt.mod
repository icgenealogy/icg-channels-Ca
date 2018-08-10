NEURON
{
  SUFFIX it2 
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

  ah = -0.1999066834808703     (/mV) 
  bh = 14.393682763381339     (1) 
  vhh = -83.76381916430927     (mV) 
  Ah = 427.79618678635643     (/ms) 
  b1h = 0.1714629516964408     (/mV) 
  c1h = 0.013460916605094264     (/mV2) 
  d1h = 0.00040029084945751115     (/mV3) 
  b2h = 0.15465040433256907     (/mV) 
  c2h = -0.0017077793160277601     (/mV2) 
  d2h = 5.444264438642098e-06     (/mV3) 

  am = 0.1351341750532069     (/mV) 
  bm = -6.621587800601443     (1) 
  vhm = -45.661714772123645     (mV) 
  Am = 2.5210620323773285     (/ms) 
  b1m = -6.976560010623724e-05     (/mV) 
  c1m = -0.00039917285282357494     (/mV2) 
  d1m = 2.4221740177371033e-06     (/mV3) 
  b2m = 0.0938898369279401     (/mV) 
  c2m = -0.001834638762245825     (/mV2) 
  d2m = 8.516655656597995e-06     (/mV3) 
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