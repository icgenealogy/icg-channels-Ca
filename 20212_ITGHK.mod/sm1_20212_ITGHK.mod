NEURON
{
  SUFFIX itGHK 
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

  ah = -0.24996475000226662     (/mV) 
  bh = 20.747259887185695     (1) 
  vhh = -90.6705139113815     (mV) 
  Ah = 138.48462562989073     (/ms) 
  b1h = -0.09601828616704017     (/mV) 
  c1h = 0.0008906403879046965     (/mV2) 
  d1h = -2.520853549364529e-06     (/mV3) 
  b2h = -0.16492234057164273     (/mV) 
  c2h = -0.024836614164963584     (/mV2) 
  d2h = -0.0015512127268020436     (/mV3) 

  am = 0.1612900180519122     (/mV) 
  bm = -9.516109813412534     (1) 
  vhm = -65.7820698773436     (mV) 
  Am = 2.8292678215178912     (/ms) 
  b1m = 0.013301214263276944     (/mV) 
  c1m = -0.0004466605355215118     (/mV2) 
  d1m = 2.9600373045117572e-06     (/mV3) 
  b2m = 0.0771455958760287     (/mV) 
  c2m = -0.0005428156392875569     (/mV2) 
  d2m = 1.2434495689125874e-06     (/mV3) 
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