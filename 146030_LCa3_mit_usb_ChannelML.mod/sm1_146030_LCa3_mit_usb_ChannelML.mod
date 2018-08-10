NEURON
{
  SUFFIX LCa3_mit_usb_ChannelML 
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

  ah = -0.12357935057447221     (/mV) 
  bh = 3.83242692694581     (1) 
  vhh = -23.552851299689898     (mV) 
  Ah = 225.03053408658315     (/ms) 
  b1h = -0.1062644548657303     (/mV) 
  c1h = 0.001367410684022844     (/mV2) 
  d1h = -5.558319212096332e-06     (/mV3) 
  b2h = -0.02057596856444799     (/mV) 
  c2h = -0.0003423966771835701     (/mV2) 
  d2h = -1.915200871843147e-06     (/mV3) 

  am = 0.14476858016657754     (/mV) 
  bm = 0.45545185265972166     (1) 
  vhm = -173.94876453379106     (mV) 
  Am = 0.023874920939284887     (/ms) 
  b1m = 0.11990846502596438     (/mV) 
  c1m = -0.0013452905534690265     (/mV2) 
  d1m = 4.9775648258015585e-06     (/mV3) 
  b2m = -0.11992045583777129     (/mV) 
  c2m = 0.0009202052490347864     (/mV2) 
  d2m = -1.8597483610240417e-06     (/mV3) 
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