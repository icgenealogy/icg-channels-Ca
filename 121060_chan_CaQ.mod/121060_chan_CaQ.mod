TITLE (Q-type calcium current for MSP Neuron)

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaQ
	USEION ca READ cai,cao WRITE ica
	RANGE minf, mtau, ica
	GLOBAL pmax
}

UNITS {
	(mA)	= (milliamp)
	(mV)	= (millivolt)
	(mM)	= (milli/liter)
        FARADAY = 96489 (coul)
        R       = 8.314 (volt-coul/degC)
}

PARAMETER {
	v		(mV)
	celsius	(degC)
	cai		(mM)
	cao		(mM)
	pmax = 6e-6	(cm/s)		
}

STATE {
	m
}

ASSIGNED {
	ica		(mA/cm2)
	mtau		(ms)
	minf
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	ica = pmax*m*m*ghk(v,cai,cao,2)
	: ica = pmax*m*m*ghk(v,0.001,cao,2)
}

DERIVATIVE state {
	rates(v)
	m'= (minf-m) / mtau
	}

INITIAL {
	rates(v)
	m = minf
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { LOCAL e, w
        w = v * (.001) * z*FARADAY / (R*(celsius+273.16))
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else : denominator is small -> Taylor series
          { e = 1-w/2 }
        ghk = - (.001) * z*FARADAY * (co-ci*exp(w)) * e
}

UNITSOFF

PROCEDURE rates(v(mV)) {
:	TABLE minf FROM -120 TO 30 WITH 70
	mtau = 0.377
	minf = 1 / (1+exp((v+9)/-6.6))
}

UNITSON 
