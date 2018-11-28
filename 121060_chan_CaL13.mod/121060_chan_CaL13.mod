TITLE (L-type LVA calcium current for MSP Neuron)

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaL13
	: USEION ca READ cai,cao 
        : USEION Ca WRITE iCa VALENCE 2
	:USEION Ca READ Cai,Cao WRITE iCa VALENCE 2
	USEION ca READ cai,cao WRITE ica
        RANGE minf, mtau, hinf, htau, pmax
	GLOBAL m_vh
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
	cai  	(mM)
	cao  	(mM)
	:Cai  	(mM)
	:Cao  	(mM)
	pmax = 4.25e-7	(cm/s)
	m_vh = -33 (mV)
}

STATE {
	m h
}

ASSIGNED {
	ica		(mA/cm2)
	mtau		(ms)
	minf
	hinf
	htau		(ms)
}

BREAKPOINT { 
	SOLVE state METHOD cnexp
	ica = pmax*m*m*h*ghk(v,cai,cao,2)
	: iCa = pmax*m*m*h*ghk(v,0.001,Cao,2)
}

DERIVATIVE state {
	rates(v)
	m'= (minf-m) / mtau
	h'= (hinf-h) / htau
	}

INITIAL {
	rates(v)
	m = minf
	h = hinf
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

PROCEDURE rates(v(mV)) { LOCAL m_alpha, m_beta
	m_alpha = 0.1194*(v+8.124)/(exp((v+8.124)/9.005)-1)
	m_beta = 2.97*exp(v/31.4)
	mtau = 1/(m_alpha+m_beta)
	htau = 14.77
	minf = 1 / (1+exp((v-m_vh)/-6.7))
	hinf = 1 / (1+exp((v+13.4)/11.9))
}

UNITSON 
