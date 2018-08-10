NEURON
{
  SUFFIX tca 
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

  aa = 0.12713677280394953     (/mV) 
  ba = -4.522313362399014     (1) 
  vha = -18.960030837771605     (mV) 
  Aa = 0.37545658333353527     (/ms) 
  b1a = -0.13857952789289857     (/mV) 
  c1a = 0.0012290919818518336     (/mV2) 
  d1a = 3.347021978348613e-05     (/mV3) 
  b2a = 0.012491123102182521     (/mV) 
  c2a = 0.00039528488942474434     (/mV2) 
  d2a = -2.6630906384988127e-06     (/mV3) 

  ab = -0.16142226598514278     (/mV) 
  bb = 10.831869403096414     (1) 
  vhb = -67.5367368529477     (mV) 
  Ab = 557.7960705248827     (/ms) 
  b1b = 0.06425448966477779     (/mV) 
  c1b = 8.289161452355732e-05     (/mV2) 
  d1b = 9.166503147290909e-07     (/mV3) 
  b2b = 0.09706039442457849     (/mV) 
  c2b = 8.871649584525644e-05     (/mV2) 
  d2b = -9.931517244934974e-07     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  aInf 
  aTau 
  bInf 
  bTau 
}

STATE
{
  a
  b
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*a*a*b
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  a' = (aInf - a) / aTau 
  b' = (bInf - b) / bTau 
}

INITIAL
{
  rates(v)
  a = aInf 
  b = bInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    aInf = 1/(1 + exp(-aa*v + ba)) 
    aTau = Aa / ( exp(-(b1a*(v-vha) + c1a*(v-vha)^2 + d1a*(v-vha)^3)) + exp((b2a*(v-vha) + c2a*(v-vha)^2 + d2a*(v-vha)^3)) ) 

    bInf = 1/(1 + exp(-ab*v + bb)) 
    bTau = Ab / ( exp(-(b1b*(v-vhb) + c1b*(v-vhb)^2 + d1b*(v-vhb)^3)) + exp((b2b*(v-vhb) + c2b*(v-vhb)^2 + d2b*(v-vhb)^3)) ) 


  UNITSON
}