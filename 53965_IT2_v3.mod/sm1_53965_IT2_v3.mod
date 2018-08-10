NEURON
{
  SUFFIX it2 
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

  ah = -0.13908985165367754     (/mV) 
  bh = 11.893819184272424     (1) 
  vhh = -91.76180465978688     (mV) 
  Ah = 1587.8243301459422     (/ms) 
  b1h = -0.1847836263222034     (/mV) 
  c1h = 0.0018750210559674257     (/mV2) 
  d1h = -5.608141574920379e-06     (/mV3) 
  b2h = -0.14085512792713595     (/mV) 
  c2h = -0.019249372479371744     (/mV2) 
  d2h = -0.0013755686838834878     (/mV3) 

  am = 0.13513504215590785     (/mV) 
  bm = -7.405398608350585     (1) 
  vhm = -77.86414906721662     (mV) 
  Am = 3.21870935724171     (/ms) 
  b1m = -0.07096453611036133     (/mV) 
  c1m = 0.0006384725419942229     (/mV2) 
  d1m = -1.8049699694379362e-06     (/mV3) 
  b2m = -0.08953260877179908     (/mV) 
  c2m = -0.0026521940007462213     (/mV2) 
  d2m = -4.930306550705305e-05     (/mV3) 
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