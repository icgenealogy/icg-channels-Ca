NEURON
{
  SUFFIX cat 
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

  ah = -0.1614849831377623     (/mV) 
  bh = 10.83544749335007     (1) 
  vhh = -78.00858967513906     (mV) 
  Ah = 63.01817365666251     (/ms) 
  b1h = 0.12321826545650787     (/mV) 
  c1h = 0.003325795072070737     (/mV2) 
  d1h = 7.10631719892144e-05     (/mV3) 
  b2h = 0.060401960991320294     (/mV) 
  c2h = -0.0005845765770804857     (/mV2) 
  d2h = 1.7335078612428651e-06     (/mV3) 

  am = 0.12713779847878115     (/mV) 
  bm = -4.522308201224841     (1) 
  vhm = -16.406391267229846     (mV) 
  Am = 2.26826253029622     (/ms) 
  b1m = -0.007912655797785093     (/mV) 
  c1m = -0.00010950786697476281     (/mV2) 
  d1m = -2.3739383552960288e-08     (/mV3) 
  b2m = 0.0890021488919847     (/mV) 
  c2m = -0.0009181649540171597     (/mV2) 
  d2m = 1.4109094605125834e-06     (/mV3) 
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