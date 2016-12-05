TITLE cat
: low-threshold Ca current for TC neurons from Williams and Stuart (2000).
: M.Migliore Jan. 2002

NEURON {
	SUFFIX cat
	USEION ca READ eca,cai,cao WRITE ica
        RANGE gbar,cai
	GLOBAL minf, hinf, mtau, htau :carev
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
} 

PARAMETER {
	gbar = 0.007   	(mho/cm2)	
								
	celsius = 22 (degC)
	v 		(mV)
	vrest=-75
	a0m=0.1
	vhalfm=55
	zetam=0.1
	gmm=0.5

	a0h=0.01
	vhalfh=10
	zetah=0.1
	gmh=0.5
	
	q10=3
	cai 	= .000050 (mM)		: initial [Ca]i = 50 nM
	cao 	= 2	(mM)		: [Ca]o = 2 mM

}



ASSIGNED {
        eca (mV)
	ica 		(mA/cm2)
	minf (mV)	mtau (ms)	 	
	hinf (mV)	htau (ms)	 	
	:carev	(mV)
}
 


STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	:carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gbar*m^3*h* (v - eca)
} 

INITIAL {
	trates(v)
	m=minf  
	h=hinf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}


PROCEDURE trates(v) {  
	LOCAL qt
        qt=q10^((celsius-22)/10)
        minf = 1/(1 + exp(-(v-(vrest+52.9))/13.6))
	mtau = betm(v)/(qt*a0m*(1+alpm(v)))

        hinf = 1/(1 + exp((v-(vrest+4.5))/11.8))
	htau = beth(v)/(qt*a0h*(1+alph(v)))
}


FUNCTION alpm(v(mV)) {
  alpm = exp(zetam*(v-vrest-vhalfm)) 
}

FUNCTION betm(v(mV)) {
  betm = exp(zetam*gmm*(v-vrest-vhalfm)) 
}

FUNCTION alph(v(mV)) {
  alph = exp(zetah*(v-vrest-vhalfh)) 
}

FUNCTION beth(v(mV)) {
  beth = exp(zetah*gmh*(v-vrest-vhalfh)) 
}
