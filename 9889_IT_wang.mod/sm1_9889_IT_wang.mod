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

  ah = -0.15873251976544875     (/mV) 
  bh = 13.57166961821516     (1) 
  vhh = -89.44070141090742     (mV) 
  Ah = 585.9574994067776     (/ms) 
  b1h = 0.08454607921558996     (/mV) 
  c1h = 0.00217447540263259     (/mV2) 
  d1h = 7.935792591553444e-06     (/mV3) 
  b2h = 0.05477390358968515     (/mV) 
  c2h = -0.0007253417231853068     (/mV2) 
  d2h = 8.789732925902869e-06     (/mV3) 

  am = 0.12820494034537952     (/mV) 
  bm = -8.333321133195877     (1) 
  vhm = -50.91271146034344     (mV) 
  Am = 1.3273708684096397     (/ms) 
  b1m = -0.07543928859165681     (/mV) 
  c1m = 0.001142605884106086     (/mV2) 
  d1m = -4.6287335747245125e-06     (/mV3) 
  b2m = -0.008612396937709536     (/mV) 
  c2m = 0.00046210863173831593     (/mV2) 
  d2m = -2.372388642196833e-06     (/mV3) 

  ad = 0.13442732294808793     (/mV) 
  bd = -10.136433178694501     (1) 
  vhd = -106.77643172888531     (mV) 
  Ad = 96.13676729137326     (/ms) 
  b1d = -0.01712424282188009     (/mV) 
  c1d = 1.3622013006855532e-05     (/mV2) 
  d1d = -3.1153312552258893e-07     (/mV3) 
  b2d = -0.09287669611086563     (/mV) 
  c2d = 0.005242304954701644     (/mV2) 
  d2d = -0.000250710873926005     (/mV3) 
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