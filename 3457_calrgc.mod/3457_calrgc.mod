TITLE CaL
: L-type Ca current for Retinal Ganglion Cell from Benison et al (2001)
: M.Migliore Nov. 2001

NEURON {
	SUFFIX calrgc
	USEION ca READ eca WRITE ica
	RANGE  gbar
	GLOBAL minf, mtau
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	
								
	eca		(mV)            : must be explicitly def. in hoc
	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ica 		(mA/cm2)
	minf		mtau (ms)
}
 

STATE { m }

BREAKPOINT {
        SOLVE states METHOD cnexp
	ica = gbar*m*m*(v - eca)
} 

INITIAL {
	trates(v)
	m=minf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
}

PROCEDURE trates(vm) {  
        LOCAL  a, b

	a = trap0(vm,3,0.061,12.5)
	b = 0.058*exp(-(vm-10)/15)
	minf = a/(a+b)
	mtau = 1/(a+b)

}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	