NEURON
{
  SUFFIX lca 
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

  ae = 0.1810995930525827     (/mV) 
  be = -0.25117903711662826     (1) 
  vhe = 3.043840158274908e-05     (mV) 
  Ae = 0.11919481258238615     (/ms) 
  b1e = -0.17895456738614104     (/mV) 
  c1e = 0.008546346996858143     (/mV2) 
  d1e = -0.0001818840354029023     (/mV3) 
  b2e = -0.17895423233881574     (/mV) 
  c2e = -0.008546372928156229     (/mV2) 
  d2e = -0.00018188569171700226     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  eInf 
  eTau 
}

STATE
{
  e
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*e*e
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  e' = (eInf - e) / eTau 
}

INITIAL
{
  rates(v)
  e = eInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    eInf = 1/(1 + exp(-ae*v + be)) 
    eTau = Ae / ( exp(-(b1e*(v-vhe) + c1e*(v-vhe)^2 + d1e*(v-vhe)^3)) + exp((b2e*(v-vhe) + c2e*(v-vhe)^2 + d2e*(v-vhe)^3)) ) 


  UNITSON
}