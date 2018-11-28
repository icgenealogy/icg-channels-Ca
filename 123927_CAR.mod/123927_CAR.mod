
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
	F = 96485 (coul)
}

PARAMETER {
	v (mV)
	celsius 		(degC)
	PcaRbar = .000044 (cm/s)
	ki=.001 (mM)
	cai=5.e-5 (mM)
	cao = 10  (mM)
	
	q10m=11.45
	q10Ampl=2.1
	q10h=3
}


NEURON {
	SUFFIX car
	USEION ca READ cai,cao WRITE ica
        RANGE PcaRbar   
        GLOBAL hinf,minf,taum,tauh
}

STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        PcaR  (cm/s) 
        minf
        hinf
        taum
        tauh
}

INITIAL {
        rates(v)
        m = minf
        h = hinf
}

UNITSOFF
BREAKPOINT {
        LOCAL qAmpl
	qAmpl = q10Ampl^((celsius - 21)/10)
	SOLVE states METHOD cnexp
	PcaR = PcaRbar*m*h
	
	ica = PcaR*qAmpl*ghk(v,cai,cao)

}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL a

        a=2*F*v/(R*(celsius+273.15)*1000)
	
        ghk=2*F/1000*(co - ci*exp(a))*func(a)
}


FUNCTION func(a) {
	if (fabs(a) < 1e-4) {
		func = -1 + a/2
	}else{
		func = a/(1-exp(a))
	}
}



DERIVATIVE states {     
        rates(v)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL alpham, f1,f2,f3,qm,qh
	
        TABLE taum, tauh, minf, hinf FROM -150 TO 150 WITH 3000
        
	qm = q10m^((celsius - 21)/10)
	qh = q10h^((celsius-21)/10)
	

	minf = 1/(1+exp(-(v+15)/5.8))
	hinf = 1/(1+exp((v+78.7)/14.5))
	
	f1=1/(1+exp(-(v+15.2)/4.29))+0.0222
	f2=15.244/(1+exp((v+13.44)/8.61))+0.511
	f3=f1*f2
        taum = f3/qm
	
        f1=1/(1+exp(-(v+49.8)/2.64))
	f2=45.11/(1+exp(v/8.92))
	f3=f1*f2+22.7
	
        tauh = f3/qh
}


UNITSON










