NEURON
{
  SUFFIX Ca 
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

  ahCa = -0.11106761854029422     (/mV) 
  bhCa = -1.111251913477313     (1) 
  vhhCa = 10.000916982055493     (mV) 
  AhCa = 999.9991043911485     (/ms) 
  b1hCa = -0.05555462625780813     (/mV) 
  c1hCa = -1.2308062226307916e-07     (/mV2) 
  d1hCa = 2.2225384518184912e-09     (/mV3) 
  b2hCa = -0.055548641739208086     (/mV) 
  c2hCa = 2.0053870047036175e-07     (/mV2) 
  d2hCa = 2.3839248668380105e-09     (/mV3) 

  amCa = 0.1666684013186313     (/mV) 
  bmCa = -3.333394285617656     (1) 
  vhmCa = 76.79676964045673     (mV) 
  AmCa = 80.77446074265109     (/ms) 
  b1mCa = -0.621661637660047     (/mV) 
  c1mCa = -0.012325587836204126     (/mV2) 
  d1mCa = -5.949966665513702e-05     (/mV3) 
  b2mCa = -0.5650946830047974     (/mV) 
  c2mCa = -0.0108538900523931     (/mV2) 
  d2mCa = -5.3454710194322954e-05     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  hCaInf 
  hCaTau 
  mCaInf 
  mCaTau 
}

STATE
{
  hCa
  mCa
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*hCa*mCa
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  hCa' = (hCaInf - hCa) / hCaTau 
  mCa' = (mCaInf - mCa) / mCaTau 
}

INITIAL
{
  rates(v)
  hCa = hCaInf 
  mCa = mCaInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    hCaInf = 1/(1 + exp(-ahCa*v + bhCa)) 
    hCaTau = AhCa / ( exp(-(b1hCa*(v-vhhCa) + c1hCa*(v-vhhCa)^2 + d1hCa*(v-vhhCa)^3)) + exp((b2hCa*(v-vhhCa) + c2hCa*(v-vhhCa)^2 + d2hCa*(v-vhhCa)^3)) ) 

    mCaInf = 1/(1 + exp(-amCa*v + bmCa)) 
    mCaTau = AmCa / ( exp(-(b1mCa*(v-vhmCa) + c1mCa*(v-vhmCa)^2 + d1mCa*(v-vhmCa)^3)) + exp((b2mCa*(v-vhmCa) + c2mCa*(v-vhmCa)^2 + d2mCa*(v-vhmCa)^3)) ) 


  UNITSON
}