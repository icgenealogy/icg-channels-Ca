TITLE (N-type calcium current for MSP Neuron)

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaN
	USEION ca READ cai,cao WRITE ica
	RANGE minf, mtau, hinf, htau
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
	pmax = 1e-5	(cm/s)		
}

CONSTANT {
	a = 0.21 (1)
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
	SOLVE state METHOD cnexp : see http://www.neuron.yale.edu/phpBB/viewtopic.php?f=28&t=592
	ica = pmax*m*m*(a*h+(1-a))*ghk(v,cai,cao,2)
	: ica = pmax*m*m*(a*h+(1-a))*ghk(v,0.001,cao,2)
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
	m_alpha = 0.1157*(v+17.19)/(exp((v+17.19)/15.22)-1)
	m_beta = 1.15*exp(v/23.82)
	mtau = 1/(m_alpha+m_beta)
	htau = 23.33
	minf = 1 / (1+exp((v+8.7)/-7.4))
	hinf = 1 / (1+exp((v+74.8)/6.5))
}

UNITSON 
