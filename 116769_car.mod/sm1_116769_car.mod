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

  ah = -0.08470638036911504     (/mV) 
  bh = 5.50631706306797     (1) 
  vhh = -102.96252914789078     (mV) 
  Ah = 395.41660113905004     (/ms) 
  b1h = 0.003181997956889208     (/mV) 
  c1h = -2.34971546165021e-05     (/mV2) 
  d1h = 3.531591063325675e-07     (/mV3) 
  b2h = 3.560156133131048e-06     (/mV) 
  c2h = 4.157959880251969e-05     (/mV2) 
  d2h = -1.2801247333727056e-07     (/mV3) 

  am = 0.14925480959984147     (/mV) 
  bm = -2.0895743072853503     (1) 
  vhm = -100.46989308302302     (mV) 
  Am = 7.219173225683626     (/ms) 
  b1m = 0.0016610687991939306     (/mV) 
  c1m = 4.422554651077385e-06     (/mV2) 
  d1m = 1.3780303446721875e-08     (/mV3) 
  b2m = 0.001624267737332303     (/mV) 
  c2m = 2.7829529139211653e-06     (/mV2) 
  d2m = -7.078707750975674e-09     (/mV3) 
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
  g = gbar*h*h*h*h*m*m*m*m*m*m*m*m*m*m*m*m
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