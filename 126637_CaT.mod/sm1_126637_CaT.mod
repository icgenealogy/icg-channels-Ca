NEURON
{
  SUFFIX CaT 
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

  ah = -0.09648604326523776     (/mV) 
  bh = 8.991203328767778     (1) 
  vhh = -93.22934940388632     (mV) 
  Ah = 400.78472858432866     (/ms) 
  b1h = -0.122608446970848     (/mV) 
  c1h = 0.0010526100654315593     (/mV2) 
  d1h = -2.811342879525347e-06     (/mV3) 
  b2h = -0.022629723987736133     (/mV) 
  c2h = -0.002239703411733503     (/mV2) 
  d2h = -0.00013264191311817918     (/mV3) 

  am = 0.1818663251874065     (/mV) 
  bm = -8.225045087109084     (1) 
  vhm = -36.679170534370215     (mV) 
  Am = 5.037495816374612     (/ms) 
  b1m = 0.002234035333377751     (/mV) 
  c1m = 0.00018938660819325     (/mV2) 
  d1m = 2.0770896507919817e-06     (/mV3) 
  b2m = 0.12112163964895713     (/mV) 
  c2m = -0.0016252762202538288     (/mV2) 
  d2m = 6.561090677207522e-06     (/mV3) 
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
  g = gbar*h*m
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