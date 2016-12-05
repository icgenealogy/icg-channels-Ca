TITLE Motoneuron N-type Calcium channels
:
: The parameters for this current come from Booth et al. JNP, 1997
:


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX N_Ca
	USEION ca READ eca WRITE ica
	RANGE gcabar
	RANGE m_inf, h_inf
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
:	gcabar  = .02 (mho/cm2)       :This is the original value
	gcabar  = .0005 (mho/cm2)
	eca	= 80	(mV)
	dt		(ms)
	tau_m	= 4	(ms)
	tau_h	= 40	(ms)
	v		(mV)
	theta_m = -20	(mV)
	kappa_m = -5  	(mV)
	theta_h = -45  	(mV)
	kappa_h = 5  	(mV)
}

STATE {
	m h
}

ASSIGNED {
	ica		(mA/cm2)
	m_inf
	h_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar * m * m * h * (v - eca)
}

DERIVATIVE states {
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) {

	m_inf = 1 / (1 + (Exp((v - theta_m)/ kappa_m)))
	h_inf = 1 / (1 + (Exp((v - theta_h)/ kappa_h)))
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(Exp(x/y)-1)
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
} 

