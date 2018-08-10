NEURON
{
  SUFFIX ca 
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

  ah = -0.10869675838110728     (/mV) 
  bh = 4.239179575942294     (1) 
  vhh = -100.0     (mV) 
  Ah = 0.17122078114114483     (/ms) 
  b1h = 5380.256888271989     (/mV) 
  c1h = 0.1958176015349293     (/mV2) 
  d1h = -1.5937848322926862e-05     (/mV3) 
  b2h = 5381.975965110726     (/mV) 
  c2h = 0.17917234094765383     (/mV2) 
  d2h = 3.150945897545833e-05     (/mV3) 

  am = 0.12047980552619031     (/mV) 
  bm = 0.36131744314752484     (1) 
  vhm = -100.0     (mV) 
  Am = 2.28     (/ms) 
  b1m = -1047.0422085239375     (/mV) 
  c1m = -0.012561218305143831     (/mV2) 
  d1m = -2.5674349241795736e-05     (/mV3) 
  b2m = -1045.208969579239     (/mV) 
  c2m = -0.03051287767709815     (/mV2) 
  d2m = 2.585379939017448e-05     (/mV3) 
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