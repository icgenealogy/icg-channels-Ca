NEURON
{
  SUFFIX mycat 
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

  ah = -0.18716955769464083     (/mV) 
  bh = 13.60811971317895     (1) 
  vhh = -39.086216874649715     (mV) 
  Ah = 119.99957722707936     (/ms) 
  b1h = -0.01015512472980993     (/mV) 
  c1h = 5.2114670454763664e-05     (/mV2) 
  d1h = -9.375236184974872e-08     (/mV3) 
  b2h = -0.01014713905532169     (/mV) 
  c2h = -5.1194078938191377e-05     (/mV2) 
  d2h = -8.742360596532905e-08     (/mV3) 

  am = 0.18715239836787015     (/mV) 
  bm = -8.79728212296968     (1) 
  vhm = -2.536103912910729     (mV) 
  Am = 40.00337433673695     (/ms) 
  b1m = -0.004669281743921934     (/mV) 
  c1m = 9.569344969303565e-06     (/mV2) 
  d1m = -1.5175633185296643e-09     (/mV3) 
  b2m = -0.004672584711603944     (/mV) 
  c2m = -1.2206850832955761e-05     (/mV2) 
  d2m = -1.3381477550033187e-08     (/mV3) 
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