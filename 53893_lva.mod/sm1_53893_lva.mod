NEURON
{
  SUFFIX lva 
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

  ah = -0.15874186002891832     (/mV) 
  bh = 13.25511197559101     (1) 
  vhh = -86.8021309221402     (mV) 
  Ah = 2869.603588121358     (/ms) 
  b1h = -0.06093996233232768     (/mV) 
  c1h = 0.0008114957050604473     (/mV2) 
  d1h = -9.758155074862721e-06     (/mV3) 
  b2h = -0.08306680308821715     (/mV) 
  c2h = -0.002202204474677861     (/mV2) 
  d2h = -1.4109465122500905e-05     (/mV3) 

  am = 0.12820497886688914     (/mV) 
  bm = -8.076913512888865     (1) 
  vhm = -67.9292297207141     (mV) 
  Am = 13.848431374219155     (/ms) 
  b1m = 0.09479718403380614     (/mV) 
  c1m = 0.0018477175561648252     (/mV2) 
  d1m = 2.9553821877878928e-05     (/mV3) 
  b2m = 0.05988754397969976     (/mV) 
  c2m = -0.0005430327933759765     (/mV2) 
  d2m = 1.5687519392669176e-06     (/mV3) 

  ad = 0.1346978563631663     (/mV) 
  bd = -9.88725920531413     (1) 
  vhd = -38.30116596664719     (mV) 
  Ad = 247.4376739602507     (/ms) 
  b1d = -0.04815255851130182     (/mV) 
  c1d = 0.0003019331372495304     (/mV2) 
  d1d = -1.544581528733546e-06     (/mV3) 
  b2d = -0.010730643359143604     (/mV) 
  c2d = -0.000622640382487268     (/mV2) 
  d2d = -7.712919232282127e-06     (/mV3) 
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
  g = gbar*h*m
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