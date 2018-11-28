: CaP high-threshold, non-inactivating Ca

NEURON {
     SUFFIX cap
     USEION ca READ eca WRITE ica
     RANGE gca, ica
}

UNITS {
     (S)  = (siemens)
     (mV) = (millivolt)
     (mA) = (milliamp)
}

PARAMETER { gca = 6e-4 (S/cm2) }

ASSIGNED {
     v       (mV)
     eca      (mV)
     ica      (mA/cm2)
}


BREAKPOINT {
     ica = gca * sinf(v) * (v - eca)
}


FUNCTION sinf (Vm (mV)) {
     UNITSOFF
     sinf = 1 / (1 + exp(-(Vm+22)/4.53)) 
     UNITSON
}



