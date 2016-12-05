NEURON {
    SUFFIX catRT03 }
NEURON {
    USEION ca READ eca WRITE ica }
ASSIGNED {
    ica
    eca (mV)
}

PARAMETER {
  :erev 		= 125.    (mV)
  gmax 		= 0.4  (S/cm2)
  vrest           = 0    (mV)

  mvhalf 	= 56.
  mkconst 	= -6.2
  exptemp 	= 37.
  mq10		= 1
  mexp 		= 2

  hvhalf 	= 80.
  hkconst 	= 4.
  hq10		= 1
  hexp 		= 1
} : end PARAMETER

INCLUDE "custom_code/inc_files/64229_boltz_cvode.inc"

FUNCTION settau(j,v) {
  if (j==0) { : m
    settau = 0.204 + .333/(exp((v+15.8)/18.2)+exp((-v-131.)/16.7))
  } else {
    if (v<-81.0) { 
      settau = 0.333*exp((v+466.)/66.6) 
    } else {
      settau = 9.32+0.333*exp((-v-21.)/10.5)
    }
  }
}

PROCEDURE iassign () { i=g*(v-eca) ica=i }
:** calRT03  -- Traub L channel
