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

  ah = -0.18189152011667437     (/mV) 
  bh = 5.839760126421066     (1) 
  vhh = -15.431014432753821     (mV) 
  Ah = 46.333188666279156     (/ms) 
  b1h = -0.0120488284935208     (/mV) 
  c1h = -1.602024335870633e-06     (/mV2) 
  d1h = 4.054063651781334e-07     (/mV3) 
  b2h = 0.024630587500232354     (/mV) 
  c2h = -0.0005565853085240663     (/mV2) 
  d2h = 2.8602289738777727e-06     (/mV3) 

  am = 0.13892087900561167     (/mV) 
  bm = -3.765709837996526     (1) 
  vhm = -26.251805460468233     (mV) 
  Am = 5.696235117359993     (/ms) 
  b1m = -0.022228646835052195     (/mV) 
  c1m = -4.890457719830979e-05     (/mV2) 
  d1m = 5.177829996992715e-07     (/mV3) 
  b2m = 0.05176840539267064     (/mV) 
  c2m = -0.000313920467048065     (/mV2) 
  d2m = -5.858544588337997e-07     (/mV3) 
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