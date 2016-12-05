TITLE N-type calcium channel 
: used in somatic and dendritic regions 
: After Borg 
:  Updated by Maria Markaki  03/12/03

NEURON {
	SUFFIX can
	USEION ca READ cai, eca WRITE ica 
 	USEION can WRITE ican VALENCE 2: to use it in can-specific pump
        RANGE gcalbar, ica, po, ican
	GLOBAL hinf, minf, s_inf
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
:	gcalbar = 0.2e-7   (cm/s)  : initialized conductance
	gcalbar = 1.0   (mho/cm2)  : initialized conductance
  	ki     = 0.025  (mM)            :test middle point of inactivation fct
  :	ki     = 0.01  (mM)            :test middle point of inactivation fct
	zetam = -3.4
	zetah = 2
	vhalfm =-21 (mV)
	vhalfh =-40 (mV)
	tm0=1.5(ms)
	th0=75(ms)
:	taumin  = 10    (ms)            : minimal value of the time cst
	taumin  = 2    (ms)            : minimal value of the time cst
}



ASSIGNED {     : parameters needed to solve DE
	v            (mV)
	celsius      (degC)
	ica          (mA/cm2)
	ican		(mA/cm2)
	po
	cai          (mM)       :5e-5 initial internal Ca++ concentration
	eca             (mV)
        minf
        hinf
	s_inf
}


FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}



STATE {	
	m 
	h 
	s
}  

INITIAL {
:	rates(v)
	rates(v,cai)
        m = minf
        h = hinf
	s = s_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	po = m*m*h
:	po = s*m*m
:	ica = gcalbar *po*h2(cai) * ghk(v,cai,cao)
 	ica = gcalbar *po*h2(cai) * (v - eca)
	ican = ica

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
:	rates(v)
	rates(v,cai)
	m' = (minf -m)/tm0
	h'=  (hinf - h)/th0
	s' = (s_inf-s)/taumin
}


:PROCEDURE rates(v (mV)) { 
PROCEDURE rates(v (mV), cai(mM)) { 
        LOCAL a, b, alpha2
        
	a = alpm(v)
	minf = 1/(1+a)
        
        b = alph(v)
	hinf = 1/(1+b)
	alpha2 = (ki/cai)^2
	s_inf = alpha2 / (alpha2 + 1)
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
