:n-type-like voltage gated calcium channel
:Migliore file Modify by Maciej Lazarewicz (mailto:mlazarew@gmu.edu) May/16/2001

TITLE n-calcium channel
: n-type calcium channel

NEURON {
	SUFFIX CAn
	USEION ca READ cai,cao WRITE ica
        RANGE gbar,ica, vshiftm, vshifth, timefactor_m, timefactor_h
        GLOBAL hinf,minf,taum,tauh, ki
    :NONSPECIFIC_CURRENT icont                                                       
}

UNITS {
	(mA) 	= 	(milliamp)
	(mV) 	= 	(millivolt)
	FARADAY =  	(faraday)  (kilocoulombs)
	R 	= 	(k-mole) (joule/degC)
	KTOMV 	= .0853 (mV/degC)
}

PARAMETER {
	v (mV)
	celsius	(degC)
    vshiftm = 0	(mV)		: voltage shift
    vshifth = 0	(mV)		: voltage shift
    timefactor_h = 1
    timefactor_m = 1
	gbar	= .0003 (mho/cm2)
	ki	= .001 	(mM)
	cai	= 5.e-5 (mM)
	cao 	= 10  	(mM)
}

STATE {	m h }

ASSIGNED {
	ica 		(mA/cm2)
    icont   (mA/cm2)
        gcan  		(mho/cm2) 
        minf
        hinf
        taum
        tauh
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcan = gbar*m*m*h*h2(cai)
	ica  = gcan*ghk(v,cai,cao)
    :icont = -ica
}

INITIAL {
	rates(v, vshiftm, vshifth)
	m = minf
	h = hinf
}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alph(v(mV)) {
	TABLE FROM -150 TO 1000 WITH 1600
	alph = 1.6e-4*exp(-v/48.4)
}

FUNCTION beth(v(mV)) {
        TABLE FROM -150 TO 1000 WITH 1600
	beth = 1/(exp((-v+39.0)/10.)+1.)
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 1000 WITH 1600
	alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
	TABLE FROM -150 TO 1000 WITH 1600
	betm = 0.046*exp(-v/20.73)
}

UNITSON

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v, vshiftm, vshifth)
        m' = (minf - m)/(taum*timefactor_m)
        h' = (hinf - h)/(tauh*timefactor_h)
}

PROCEDURE rates(v (mV), vshiftm (mV), vshifth (mV)) { :callable from hoc
        LOCAL a

        a    = alpm(v-vshiftm)
        taum = 1/(a + betm(v-vshiftm))
        minf = a*taum

        a    = alph(v-vshifth)
        tauh = 1/(a + beth(v-vshifth))
        hinf = a*tauh
}










