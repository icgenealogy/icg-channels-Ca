NEURON
{
  SUFFIX calrgc 
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

  am = 0.11586119982567059     (/mV) 
  bm = -1.8195531188848688     (1) 
  vhm = -13.986269969896401     (mV) 
  Am = 3.1134198427494417     (/ms) 
  b1m = 0.05942907846012127     (/mV) 
  c1m = -0.00016994313837112243     (/mV2) 
  d1m = -1.8014213169472198e-06     (/mV3) 
  b2m = 0.05462066027672934     (/mV) 
  c2m = -0.0004019675893723442     (/mV2) 
  d2m = 1.2903359344717288e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  mInf 
  mTau 
}

STATE
{
  m
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*m*m
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  m' = (mInf - m) / mTau 
}

INITIAL
{
  rates(v)
  m = mInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    mInf = 1/(1 + exp(-am*v + bm)) 
    mTau = Am / ( exp(-(b1m*(v-vhm) + c1m*(v-vhm)^2 + d1m*(v-vhm)^3)) + exp((b2m*(v-vhm) + c2m*(v-vhm)^2 + d2m*(v-vhm)^3)) ) 


  UNITSON
}