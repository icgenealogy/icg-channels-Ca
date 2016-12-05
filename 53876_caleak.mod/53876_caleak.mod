COMMENT
This file, caleak.mod, implements the Ca leak current G_Ca(leak) in 
Quadroni and Knopfel 1994 table 2
ENDCOMMENT

NEURON {
  SUFFIX caleak
  : NONSPECIFIC_CURRENT i
  USEION ca WRITE  ica
  RANGE i, Erev, g
}

PARAMETER {
  g = 1.33e-6 (siemens/cm2) < 0, 1e9 >
  Erev = 80 (millivolt)
}

ASSIGNED {
  i (milliamp/cm2)
  ica (milliamp/cm2)
  v (millivolt)
}

BREAKPOINT { 
  ica = g * (v - Erev)
  i = ica	: i is used for diagnostic purposes only (visible in hoc)
}
