TITLE T-type calcium channel

COMMENT
T-Type Ca2+ channel
From: Migliore et al, 1995; based on Jaffe et al, 1994
Updates:
20100910-MJCASE-documented
ENDCOMMENT

VERBATIM
#include <stdlib.h> /* 	Include this library so that the following
						(innocuous) warning does not appear:
						 In function '_thread_cleanup':
						 warning: incompatible implicit declaration of 
						          built-in function 'free'  */
ENDVERBATIM

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER {
	v (mV) 					: membrane potential
      celsius (degC) : temperature - set in hoc; default is 6.3
	gmax=.003 (mho/cm2)	: conductance flux (bar(?))
	cai (mM)				: intracellular Ca2+ concentration
	cao (mM)				: extracellular Ca2+ concentration
}

NEURON {
	SUFFIX ch_CavT				: The name of the mechanism
	USEION ca READ cai, cao WRITE ica VALENCE 2
	RANGE gmax, cai, ica
	RANGE myi
	THREADSAFE
}

STATE {
	m h 	: m = , h = 
}

ASSIGNED {			: assigned (where?)
	ica (mA/cm2)	: current flux
    g (mho/cm2)	: conductance flux
	:etca (mV)		: reversal potential
	myi (mA/cm2)
}

: verbatim blocks are not thread safe (perhaps related, this mechanism cannot be used with cvode)
INITIAL {
      m = minf(v)
      h = hinf(v)
	VERBATIM			
	cai=_ion_cai;
	ENDVERBATIM
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*m*m*h
	ica = g*ghk(v,cai,cao)
	myi = ica

}

DERIVATIVE states {	: exact when v held constant
	m' = (minf(v) - m)/m_tau(v)
	h' = (hinf(v) - h)/h_tau(v)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION hinf(v(mV)) {
	LOCAL a,b
	TABLE FROM -150 TO 150 WITH 200
	a = 1.e-6*exp(-v/16.26)
	b = 1/(exp((-v+29.79)/10)+1)
	hinf = a/(a+b)
}

FUNCTION minf(v(mV)) {
	LOCAL a,b
	TABLE FROM -150 TO 150 WITH 200
        
	a = 0.2*(-1.0*v+19.26)/(exp((-1.0*v+19.26)/10.0)-1.0)
	b = 0.009*exp(-v/22.03)
	minf = a/(a+b)
}

FUNCTION m_tau(v(mV)) (ms) {
	LOCAL a,b
	TABLE FROM -150 TO 150 WITH 200
	a = 0.2*(-1.0*v+19.26)/(exp((-1.0*v+19.26)/10.0)-1.0)
	b = 0.009*exp(-v/22.03)
	m_tau = 1/(a+b)
}

FUNCTION h_tau(v(mV)) (ms) {
	LOCAL a,b
        TABLE FROM -150 TO 150 WITH 200
	a = 1.e-6*exp(-v/16.26)
	b = 1/(exp((-v+29.79)/10.)+1.)
	h_tau = 1/(a+b)
}
