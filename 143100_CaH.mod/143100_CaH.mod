TITLE high-voltage-activated calcium channel for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (mM) = (milli/liter)
}

NEURON {
 SUFFIX CaH
 USEION ca READ eca WRITE ica
 RANGE gmax, iCaH
}

PARAMETER {
 v (mV)
 dt (ms)
 :cai (mM)
 :cao (mM)
 gmax  = 0.001 (mho/cm2)
 iCaH  = 0.0 (mA/cm2)
 :e = 130 (mV)

 theta_m = -20.0 (mV)
 k_m = 7.0 (mV)
 taum = 0.2 (ms)
}

STATE {
 m
}

ASSIGNED { 
 eca (mV)
 ica (mA/cm2)
 minf
}

BREAKPOINT {
 SOLVE states METHOD cnexp
 ica  = gmax*m*(v-eca)
 iCaH = ica
}

UNITSOFF

INITIAL {
 settables(v)
 m = minf
}

DERIVATIVE states {  
 settables(v)
 m' = (minf - m)/taum
}

PROCEDURE settables(v) {
        TABLE minf FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
}

UNITSON
