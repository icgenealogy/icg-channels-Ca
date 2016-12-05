TITLE Motoneuron HVA L-type Calcium channels
:
: The parameters for this current come from V. Booth et al. J Neurophysiol 78:3371-3385, 1997
: Iterative equations
:	I shifted the activation of this current to the same range as N_Ca.mod 
:



NEURON {
	SUFFIX L_HVA_Ca
	USEION ca READ eca WRITE ica
	RANGE gcabar
	RANGE m_inf
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gcabar  = 0.0003 (mho/cm2)  :original = 7e-05
	eca	= 80	(mV)
:	celcius = 36	(degC)
	dt		(ms)
	tau_m	= 20	(ms)
	v		(mV)
	theta_m = -10	(mV)
	kappa_m = -6    (mV)	:original = -7
}

STATE {
	m 
}

ASSIGNED {
	ica		(mA/cm2)
	m_inf
	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar * m * (v - eca)  :I have tried this as m*m also
}

DERIVATIVE states {
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {

:
:  Q10 was assumed to be 3 for both currents
:
:	tadj = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV)) {

	m_inf = 1 / (1 + (Exp((v - theta_m)/ kappa_m))): / tadj

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
