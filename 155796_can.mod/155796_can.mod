TITLE N-type Ca current 
: Taken from Retinal Ganglion Cell from Benison et al (2001), M.Migliore Nov. 2001

NEURON {
	SUFFIX can
	USEION ca READ eca WRITE ica
	RANGE  gbar
	GLOBAL minf, hinf, mtau, htau
	THREADSAFE minf, hinf, mtau, htau
}

PARAMETER {
    gbar = 4e-2   	(mho/cm2)
    mafac = 0.1
    mamid = -15
    maslope = 5
    mbfac = 0.5
    mbmid = -5
    mbslope = 10
    
    hafac = 0.01
    hamid = 50
    haslope = 15
    hbfac = 0.2
    hbmid = 40
    hbslope = 17
    
    eca		(mV)            
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
	minf 		hinf 		
	mtau (ms)	htau (ms) 	
}
 

STATE { m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ica = gbar*m*m*h * (v - eca)
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

PROCEDURE trates(vm) {  
        LOCAL  a, b
	
	:a = trap0(vm,0,0.1,10)
	a = trap0(vm,mamid,mafac,maslope)
	b = mbfac*exp(-(vm+mbmid)/mbslope)
	minf = a/(a+b)
	mtau = 1/(a+b)

	a = hafac*exp(-(vm+hamid)/haslope)
	b = hbfac/(1+exp(-(vm+hbmid)/hbslope))
	hinf = a/(a+b)
	htau = 1/(a+b)
}

FUNCTION trap0(v,th,a,q) {
	if (fabs(v-th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}

