NEURON
{
  SUFFIX CAn 
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

  ah = -0.12011681526223732     (/mV) 
  bh = 4.822370831391353     (1) 
  vhh = -43.03782240039915     (mV) 
  Ah = 3016.3638026064623     (/ms) 
  b1h = -0.08600827128713115     (/mV) 
  c1h = -0.00039637589201494625     (/mV2) 
  d1h = 4.314731968248707e-06     (/mV3) 
  b2h = -0.02779829945846945     (/mV) 
  c2h = -0.00012367181601081478     (/mV2) 
  d2h = -7.615960836517164e-07     (/mV3) 

  am = 0.12504121178315392     (/mV) 
  bm = -2.6876326194330935     (1) 
  vhm = -18.218308716323236     (mV) 
  Am = 7.143044248406778     (/ms) 
  b1m = 0.03425999436921886     (/mV) 
  c1m = -0.0002769662100279975     (/mV2) 
  d1m = -1.846807889211035e-06     (/mV3) 
  b2m = 0.0882227174590348     (/mV) 
  c2m = -0.0005678614091559835     (/mV2) 
  d2m = 1.197248500394733e-06     (/mV3) 
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
  g = gbar*h*h*m*m
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