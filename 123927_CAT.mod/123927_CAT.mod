TITLE T-calcium channel

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
	celsius = 6.3	(degC)
: wie gross ist PcaTbar
	PcaTbar = .000011 (cm/s)
	cai (mM)
	cao (mM)
	q10Ampl=3.3
	q10m=3.55
	q10h=2.8
}


NEURON {
	SUFFIX cat
	USEION ca READ cai,cao WRITE ica
        RANGE PcaTbar,cai
}

STATE {
	m h 
}

ASSIGNED {
	ica (mA/cm2)
        PcaT (cm/s)
}

INITIAL {
      m = minf(v)
      h = hinf(v)
}

UNITSOFF
BREAKPOINT {
	SOLVE states METHOD cnexp
	PcaT = PcaTbar*m*m*h
	ica = PcaT*ghk(v,cai,cao)

}

DERIVATIVE states {	
	m' = (minf(v) - m)/m_tau(v)
	h' = (hinf(v) - h)/h_tau(v)
}




FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL a, qtAmpl
	
	qtAmpl=q10Ampl^((celsius-23)/10)

        a=2*F*v/(R*(celsius+273.15)*1000)
	
        ghk=qtAmpl*2*F/1000*(co - ci*exp(a))*func(a)
}


FUNCTION func(a) {
	if (fabs(a) < 1e-4) {
		func = -1 + a/2
	}else{
		func = a/(1-exp(a))
	}
}


FUNCTION hinf(v(mV)) 
{       TABLE FROM -150 TO 150 WITH 3000 :Mitti 
	hinf = 1/(1+exp((v+72)/3.7))
}

FUNCTION minf(v(mV)) {
	TABLE FROM -150 TO 150 WITH 3000 :Mitti
        minf = (1/(1+exp(-(v+31.4)/8.8)))^0.5
}

FUNCTION m_tau(v(mV)) (ms) {
	LOCAL f1,f2, qtm
	
        TABLE FROM -150 TO 150 WITH 3000 :Mitti
        
	qtm=q10m^((celsius-23)/10)
	
	f1=1/(1+exp(-(v-7.63)/28.47))+0.01
	f2=62.82/(1+exp((v+37.02)/5.27))+3.78

	m_tau=f1*f2/qtm
}

FUNCTION h_tau(v(mV)) (ms) {
	LOCAL alphah, localhinf,qth
	
        TABLE FROM -150 TO 150 WITH 3000 :Mitti
        
	qth=q10h^((celsius-23)/10)
	
	localhinf = 1/(1+exp((v+72)/3.7))
	
	alphah=0.0021/(1+exp((v+65.77)/4.32))
	
	h_tau = localhinf/(qth*alphah)
}

UNITSON







