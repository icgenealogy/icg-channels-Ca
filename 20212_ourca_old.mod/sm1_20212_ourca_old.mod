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

  ah = -1.0003977638594672     (/mV) 
  bh = 52.022346382004805     (1) 
  vhh = -52.266501416394995     (mV) 
  Ah = 19.98632214374881     (/ms) 
  b1h = 7.307096995663666e-05     (/mV) 
  c1h = -2.593881909271725e-05     (/mV2) 
  d1h = 1.5998558211418775e-07     (/mV3) 
  b2h = 7.977538370802019e-05     (/mV) 
  c2h = -2.942848604552289e-05     (/mV2) 
  d2h = 1.816650741676165e-07     (/mV3) 

  am = 0.33323889139461504     (/mV) 
  bm = -16.32917871248497     (1) 
  vhm = -38.75944733620522     (mV) 
  Am = 99.9969779928933     (/ms) 
  b1m = -5.875910346003366e-07     (/mV) 
  c1m = -7.85633119315516e-06     (/mV2) 
  d1m = 4.903999227912013e-08     (/mV3) 
  b2m = 3.3815085438427314e-07     (/mV) 
  c2m = -8.240464129005623e-06     (/mV2) 
  d2m = 5.142590537749804e-08     (/mV3) 
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
  g = gbar*h*m*m*m
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