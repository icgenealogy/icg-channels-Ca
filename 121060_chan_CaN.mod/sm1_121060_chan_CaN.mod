NEURON
{
  SUFFIX CaN 
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

  ah = -0.15384151152085157     (/mV) 
  bh = 11.507389396795789     (1) 
  vhh = -31.159213127286694     (mV) 
  Ah = 46.6800000038824     (/ms) 
  b1h = 7.03201715367646e-10     (/mV) 
  c1h = 2.197249832290033e-10     (/mV2) 
  d1h = 6.591365126522875e-11     (/mV3) 
  b2h = 7.49912096221727e-10     (/mV) 
  c2h = 2.1988452723839914e-10     (/mV2) 
  d2h = 6.590114652087247e-11     (/mV3) 

  am = 0.135134931315481     (/mV) 
  bm = -1.1756617307085253     (1) 
  vhm = -8.722691913208475     (mV) 
  Am = 0.9564507548373041     (/ms) 
  b1m = 0.04175740422666518     (/mV) 
  c1m = 0.0002670274278823318     (/mV2) 
  d1m = 5.978923436616842e-07     (/mV3) 
  b2m = 0.02922812650148947     (/mV) 
  c2m = 0.00020427965898678138     (/mV2) 
  d2m = -1.2523426805402292e-06     (/mV3) 
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