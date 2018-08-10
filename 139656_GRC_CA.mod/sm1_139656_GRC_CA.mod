NEURON
{
  SUFFIX GRC_CA 
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

  au = -0.06696662245789706     (/mV) 
  bu = 3.2144123657668393     (1) 
  vhu = 24.78182203197156     (mV) 
  Au = 100.09553850145554     (/ms) 
  b1u = -0.007360148335931063     (/mV) 
  c1u = -0.00011049276058956579     (/mV2) 
  d1u = 1.0758300133350507e-06     (/mV3) 
  b2u = 0.007360240608086096     (/mV) 
  c2u = -0.00039140951074812723     (/mV2) 
  d2u = -1.013433902666373e-06     (/mV3) 

  as = 0.10194438846642868     (/mV) 
  bs = -2.039455621501622     (1) 
  vhs = -11.48328714201013     (mV) 
  As = 1.4765738526734085     (/ms) 
  b1s = -0.08786635937129204     (/mV) 
  c1s = 0.0006580440821973302     (/mV2) 
  d1s = 1.0209320118440365e-07     (/mV3) 
  b2s = -0.017742273708559998     (/mV) 
  c2s = 0.00025963084731378815     (/mV2) 
  d2s = 1.0427655966104453e-06     (/mV3) 
}

ASSIGNED
{
  v	(mV)
  eca	(mV)
  ica	(mA/cm2)
  g	(S/cm2)
  celsius	(degC)
  uInf 
  uTau 
  sInf 
  sTau 
}

STATE
{
  u
  s
}

BREAKPOINT
{
  SOLVE states METHOD cnexp
  g = gbar*u*s*s
  ica = g*(v-eca)
}

DERIVATIVE states
{
  rates(v)
  u' = (uInf - u) / uTau 
  s' = (sInf - s) / sTau 
}

INITIAL
{
  rates(v)
  u = uInf 
  s = sInf 
}

PROCEDURE rates(v(mV))
{
  UNITSOFF

    uInf = 1/(1 + exp(-au*v + bu)) 
    uTau = Au / ( exp(-(b1u*(v-vhu) + c1u*(v-vhu)^2 + d1u*(v-vhu)^3)) + exp((b2u*(v-vhu) + c2u*(v-vhu)^2 + d2u*(v-vhu)^3)) ) 

    sInf = 1/(1 + exp(-as*v + bs)) 
    sTau = As / ( exp(-(b1s*(v-vhs) + c1s*(v-vhs)^2 + d1s*(v-vhs)^3)) + exp((b2s*(v-vhs) + c2s*(v-vhs)^2 + d2s*(v-vhs)^3)) ) 


  UNITSON
}