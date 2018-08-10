NEURON
{
  SUFFIX cachan 
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

  afhva = -0.19999960883568896     (/mV) 
  bfhva = 9.599976013059178     (1) 
  vhfhva = -14.762096907039341     (mV) 
  Afhva = 0.9609701517885835     (/ms) 
  b1fhva = -0.0068194863559649195     (/mV) 
  c1fhva = -0.00010718026542368478     (/mV2) 
  d1fhva = 1.2118901329087776e-06     (/mV3) 
  b2fhva = 0.006819358691710315     (/mV) 
  c2fhva = -0.001225265940878558     (/mV2) 
  d2fhva = 9.689008993017756e-06     (/mV3) 

  adl = 0.19999954192003794     (/mV) 
  bdl = -8.999978756453212     (1) 
  vhdl = -38.166638567841105     (mV) 
  Adl = 7.520613031752074     (/ms) 
  b1dl = -0.047487190723764824     (/mV) 
  c1dl = -0.0003041913269413491     (/mV2) 
  d1dl = 9.81708392877752e-06     (/mV3) 
  b2dl = 0.13327353364185326     (/mV) 
  c2dl = -0.001639908478916086     (/mV2) 
  d2dl = 6.188255185457596e-06     (/mV3) 

  adhva = 0.09999987461300205     (/mV) 
  bdhva = -0.9999908360392831     (1) 
  vhdhva = -44.93822129239289     (mV) 
  Adhva = 0.1763715111219629     (/ms) 
  b1dhva = 0.0033908571590204333     (/mV) 
  c1dhva = -0.0003603841742976498     (/mV2) 
  d1dhva = 2.1608486687163117e-06     (/mV3) 
  b2dhva = 0.07671530045872495     (/mV) 
  c2dhva = -0.0019361822434065634     (/mV2) 
  d2dhva = 9.851201027466228e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  fhvaInf 
  fhvaTau 
  dlInf 
  dlTau 
  dhvaInf 
  dhvaTau 
}

STATE
{
  fhva
  dl
  dhva
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*fhva*dl*dhva
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  fhva' = (fhvaInf - fhva) / fhvaTau 
  dl' = (dlInf - dl) / dlTau 
  dhva' = (dhvaInf - dhva) / dhvaTau 
}

INITIAL
{
  rates(v)
  fhva = fhvaInf 
  dl = dlInf 
  dhva = dhvaInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    fhvaInf = 1/(1 + exp(-afhva*v + bfhva)) 
    fhvaTau = Afhva / ( exp(-(b1fhva*(v-vhfhva) + c1fhva*(v-vhfhva)^2 + d1fhva*(v-vhfhva)^3)) + exp((b2fhva*(v-vhfhva) + c2fhva*(v-vhfhva)^2 + d2fhva*(v-vhfhva)^3)) ) 

    dlInf = 1/(1 + exp(-adl*v + bdl)) 
    dlTau = Adl / ( exp(-(b1dl*(v-vhdl) + c1dl*(v-vhdl)^2 + d1dl*(v-vhdl)^3)) + exp((b2dl*(v-vhdl) + c2dl*(v-vhdl)^2 + d2dl*(v-vhdl)^3)) ) 

    dhvaInf = 1/(1 + exp(-adhva*v + bdhva)) 
    dhvaTau = Adhva / ( exp(-(b1dhva*(v-vhdhva) + c1dhva*(v-vhdhva)^2 + d1dhva*(v-vhdhva)^3)) + exp((b2dhva*(v-vhdhva) + c2dhva*(v-vhdhva)^2 + d2dhva*(v-vhdhva)^3)) ) 


  UNITSON
}