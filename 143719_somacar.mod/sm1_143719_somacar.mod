NEURON
{
  SUFFIX somacar 
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

  ah = -1.0006607274261317     (/mV) 
  bh = 62.044125848453106     (1) 
  vhh = -126.05283695746009     (mV) 
  Ah = 5.02906691835457     (/ms) 
  b1h = -7.239789706064386e-05     (/mV) 
  c1h = 4.448727375886119e-07     (/mV2) 
  d1h = -8.848853240677416e-10     (/mV3) 
  b2h = -0.4351462455999367     (/mV) 
  c2h = 0.009316420041156941     (/mV2) 
  d2h = -6.076288122405101e-05     (/mV3) 

  am = 0.33321735759066057     (/mV) 
  bm = -19.993826076280477     (1) 
  vhm = -20.772845080799954     (mV) 
  Am = 200.02868482465854     (/ms) 
  b1m = -0.0036042500991813314     (/mV) 
  c1m = -7.013057296158942e-07     (/mV2) 
  d1m = 2.2871712728121246e-08     (/mV3) 
  b2m = -0.0035887265735091026     (/mV) 
  c2m = -1.3913138301635415e-05     (/mV2) 
  d2m = -3.1785413018155005e-08     (/mV3) 
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