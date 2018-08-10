NEURON
{
  SUFFIX Ca_HVA 
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

  ah = -0.04675273845333452     (/mV) 
  bh = 2.690546822126059     (1) 
  vhh = -7.536008748395804     (mV) 
  Ah = 491.6891995687899     (/ms) 
  b1h = -0.012635583778398003     (/mV) 
  c1h = 0.00020219384351154418     (/mV2) 
  d1h = 2.129941775237684e-06     (/mV3) 
  b2h = 0.012634776676964304     (/mV) 
  c2h = 7.408877744538425e-05     (/mV2) 
  d2h = -8.647718338196585e-07     (/mV3) 

  am = 0.2189075135732005     (/mV) 
  bm = -7.232253832089235     (1) 
  vhm = -36.55026033365949     (mV) 
  Am = 13.7568025239211     (/ms) 
  b1m = 0.09332592612668503     (/mV) 
  c1m = 0.001196640466214618     (/mV2) 
  d1m = 1.2279911041764199e-05     (/mV3) 
  b2m = 0.13792877538942075     (/mV) 
  c2m = -0.0016916749459321505     (/mV2) 
  d2m = 7.1320946425071065e-06     (/mV3) 
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