TITLE t-type calcium channel with high threshold for activation
: used in somatic and dendritic regions 
:
: 
: Updated to use CVode --Carl Gold 08/10/03


NEURON {
	SUFFIX cat
	USEION ca READ cai, eca WRITE ica   
        :USEION Ca WRITE iCa VALENCE 2
        : The T-current does not activate calcium-dependent K-currents
        :RANGE gbar, iCa
        RANGE gbar, ica
	GLOBAL hinf, minf
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {           :parameters that can be entered when function is called in cell-setup 
:	gbar = 0.1e-7   (cm/s)  : initialized conductance
	gbar = 1.0   (mho/cm2)  : initialized conductance
	zetam = -3
	zetah = 5.2
	vhalfm =-36 (mV)
	vhalfh =-68 (mV)
	tm0=1.5(ms)
	th0=10(ms)
}



ASSIGNED {     : parameters needed to solve DE
	v            (mV)
	celsius      (degC)
:	iCa          (mA/cm2)
	ica          (mA/cm2)
	cai          (mM)       :5e-5 initial internal Ca++ concentration
	eca          (mV)       : initial external Ca++ concentration
        minf
        hinf
}


STATE {	
	m 
	h 
}  

INITIAL {
	rates(v)
        m = minf
        h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp

:	ecat = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
:	iCa = gbar*m*m*h*(v-eca)	: dummy calcium current induced by this channel
	ica = gbar*m*m*h*(v-eca)	: dummy calcium current induced by this channel

}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
	eci = ci*efun(-z)
	:high cao charge moves inward
	:negative potential charge moves inward
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}


DERIVATIVE states {
	rates(v)
	m' = (minf -m)/tm0
	h'=  (hinf - h)/th0
}


PROCEDURE rates(v (mV)) { 
        LOCAL a, b
        
	a = alpm(v)
	minf = 1/(1+a)
        
        b = alph(v)
	hinf = 1/(1+b)
}



FUNCTION alpm(v(mV)) {
UNITSOFF
  alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION alph(v(mV)) {
UNITSOFF
  alph = exp(1.e-3*zetah*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

