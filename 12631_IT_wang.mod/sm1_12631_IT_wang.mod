NEURON
{
  SUFFIX it 
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

  ah = -0.15873251976545766     (/mV) 
  bh = 13.571669618215944     (1) 
  vhh = -89.44747316976458     (mV) 
  Ah = 585.3000389241638     (/ms) 
  b1h = 0.08455828038627854     (/mV) 
  c1h = 0.0021759754449979145     (/mV2) 
  d1h = 8.01284055464636e-06     (/mV3) 
  b2h = 0.05474213083892965     (/mV) 
  c2h = -0.0007249953865484895     (/mV2) 
  d2h = 8.786454921927438e-06     (/mV3) 

  am = 0.12820494034537266     (/mV) 
  bm = -8.333321133195435     (1) 
  vhm = -50.91271146034344     (mV) 
  Am = 1.3273708684096397     (/ms) 
  b1m = -0.07543928859165681     (/mV) 
  c1m = 0.001142605884106086     (/mV2) 
  d1m = -4.6287335747245125e-06     (/mV3) 
  b2m = -0.008612396937709536     (/mV) 
  c2m = 0.00046210863173831593     (/mV2) 
  d2m = -2.372388642196833e-06     (/mV3) 

  ad = 0.13442730509796103     (/mV) 
  bd = -10.136431701665373     (1) 
  vhd = -106.81006143399885     (mV) 
  Ad = 96.19988059295433     (/ms) 
  b1d = -0.017128098959107923     (/mV) 
  c1d = 1.3686178751829333e-05     (/mV2) 
  d1d = -3.1164514658365153e-07     (/mV3) 
  b2d = -0.09214315475680844     (/mV) 
  c2d = 0.0051771672898317     (/mV2) 
  d2d = -0.0002485033006496047     (/mV3) 
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
  dInf 
  dTau 
}

STATE
{
  h
  m
  d
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
  d' = (dInf - d) / dTau 
}

INITIAL
{
  rates(v)
  h = hInf 
  m = mInf 
  d = dInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hInf = 1/(1 + exp(-ah*v + bh)) 
    hTau = Ah / ( exp(-(b1h*(v-vhh) + c1h*(v-vhh)^2 + d1h*(v-vhh)^3)) + exp((b2h*(v-vhh) + c2h*(v-vhh)^2 + d2h*(v-vhh)^3)) ) 

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 

    dInf = 1/(1 + exp(-ad*v + bd)) 
    dTau = Ad / ( exp(-(b1d*(v-vhd) + c1d*(v-vhd)^2 + d1d*(v-vhd)^3)) + exp((b2d*(v-vhd) + c2d*(v-vhd)^2 + d2d*(v-vhd)^3)) ) 


  UNITSON
}