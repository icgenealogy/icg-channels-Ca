TITLE L-type calcium channel for Tiger Salamander Bipolar cell
:
: Modified from Fohlmeister et al, 1990, Brain Res 510, 343-345
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX tsbp
	USEION ca READ cai, eca, cao WRITE ica
	RANGE gcabar
	RANGE c_inf
	RANGE tau_c
	RANGE c_exp

}


UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	gcabar	= 0.002	(mho/cm2)
	eca		(mV)
	cao	= 1.8	(mM)
	cai     = 0.0001 (mM)
	dt              (ms)
	v               (mV)

}

STATE {
	c 
}

INITIAL {
: The initial values were determined at a resting value of -66.3232 mV in a single-compartment
:	c = 0.0016
: at -60 mV
        c = 0.0038
}

ASSIGNED {
	ica	(mA/cm2)
	c_inf
	tau_c
	c_exp

}

BREAKPOINT {
	SOLVE states
	ica = gcabar * c*c*c * (v - eca)

}

PROCEDURE states() {	: exact when v held constant
	evaluate_fct(v)
	c = c + c_exp * (c_inf - c)

	VERBATIM
	return 0;
	ENDVERBATIM

}

UNITSOFF

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b
	
:CA channel
: a = (-0.3 * (v+10)) / ((exp(-0.1*(v+10))) - 1) : for Goldfish Bipolar
 a = (-0.3 * (v+70)) / ((exp(-0.1*(v+70))) - 1)
 b = 10 * (exp((-1*(v + 38))/9))


	tau_c = 1 / (a + b)
	c_inf = a * tau_c

: State vars to inifinity
	c_exp = 1 - exp(-dt/tau_c)

}

UNITSON
