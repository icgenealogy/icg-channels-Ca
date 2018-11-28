
COMMENT

gBoltzT.mod

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Calcium
	USEION ca READ eca WRITE ica
	RANGE n, gca, gbar, ninf, nexp
	GLOBAL v12, vSlope, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 150   	(pS/um2)	: 0.03 mho/cm2
	v 		(mV)
								
	v12 = -17.3	(mV)
	vSlope = 11.3	(mV)
	tau = 3
} 


ASSIGNED {
	ica 		(mA/cm2)
	gca		(pS/um2)
	eca		(mV)
	ninf
	nexp
}
 

STATE { n }

INITIAL { 
	states()
	
}

BREAKPOINT {
        SOLVE states
	gca = gbar*n
	ica = (1e-4) * gca * (v - eca)
} 

PROCEDURE states() {   : Computes state variable n 

	nexp = 1-exp(-dt/tau)
	ninf = 1/(1+exp(-(v-v12)/vSlope))
	n = n + nexp*(ninf-n)
}
