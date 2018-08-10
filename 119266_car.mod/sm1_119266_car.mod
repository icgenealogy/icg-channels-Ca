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

  ah = -0.07482823496617237     (/mV) 
  bh = 2.993560188998912     (1) 
  vhh = -3.3502027740126628     (mV) 
  Ah = 149.9391503251252     (/ms) 
  b1h = -0.016166029409379622     (/mV) 
  c1h = 0.0001314168288645116     (/mV2) 
  d1h = -3.7165551872531603e-07     (/mV3) 
  b2h = -0.016156805971557588     (/mV) 
  c2h = -0.00013126626820322254     (/mV2) 
  d2h = -3.7111807112733975e-07     (/mV3) 

  am = 0.12724346197338585     (/mV) 
  bm = -2.6734560505082148     (1) 
  vhm = -49.99906804498918     (mV) 
  Am = 3.0200000007215264     (/ms) 
  b1m = 2.3495020870022847e-08     (/mV) 
  c1m = -1.8506016871513104e-09     (/mV2) 
  d1m = 4.5661881409676185e-11     (/mV3) 
  b2m = 2.351113689296779e-08     (/mV) 
  c2m = -1.8501187228770067e-09     (/mV2) 
  d2m = 4.565469124883422e-11     (/mV3) 
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