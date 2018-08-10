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

  ah = -0.16109705096192442     (/mV) 
  bh = 10.812724632968864     (1) 
  vhh = -54.225789465445054     (mV) 
  Ah = 3953.2790308739154     (/ms) 
  b1h = -0.08934351582460028     (/mV) 
  c1h = -0.00033618297448365667     (/mV2) 
  d1h = 3.7842024947564564e-06     (/mV3) 
  b2h = -0.06865568150050366     (/mV) 
  c2h = -0.00018077880387852764     (/mV2) 
  d2h = -1.6028678737259995e-06     (/mV3) 

  am = 0.12700991051219213     (/mV) 
  bm = -4.517755062842831     (1) 
  vhm = -32.82112833391632     (mV) 
  Am = 20.58215483844669     (/ms) 
  b1m = -0.09149488981891066     (/mV) 
  c1m = 0.00040511002832993966     (/mV2) 
  d1m = -1.566070624370007e-07     (/mV3) 
  b2m = -0.033087950232754515     (/mV) 
  c2m = 0.0002605525631157622     (/mV2) 
  d2m = 1.839401206205263e-06     (/mV3) 
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