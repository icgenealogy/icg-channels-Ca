UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mS) = (millisiemens)
}

NEURON {
  SUFFIX icalts
  USEION ca WRITE ica
  RANGE gca,eca
}

PARAMETER {
  gca = 1    (mS/cm2)
  eca = 120  (mV)
}
    
ASSIGNED { 
  ica (mA/cm2)    
  v   (mV)
}

PROCEDURE iassign () { ica = (1e-3) * gca * mcainf(v)^2 * (v-eca) }

INITIAL {
  iassign()
}

BREAKPOINT { iassign() }

FUNCTION mcainf(v(mV)) { mcainf = fun2(v, -20, 1,  -9)*1(ms) }

INCLUDE "aux_fun.inc"
