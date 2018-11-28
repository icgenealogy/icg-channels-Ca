TITLE Motoneuron L-type Calcium channels
:
: The parameters for this current come from V. Booth et al. J Neurophysiol 78:3371-3385, 1997
: Iterative equations
: Modified by RP to provide calcium to a separate pool (caL)and to have adjustable equilibrium
: potential vca


NEURON {
	SUFFIX L_Ca
	:USEION caL READ ecaL WRITE icaL VALENCE 2
	:USEION ca READ eca
	USEION ca READ eca WRITE ica
	RANGE gcabar,ica,m_inf,m
	GLOBAL vca,theta_m,kappa_m,tau_m, eca
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gcabar  = 0.0003  (mho/cm2)
:	celcius = 36	(degC)
	dt		(ms)
	tau_m	= 20	(ms)
	v		(mV)
    vca=80		(mV)
	theta_m = -30   (mV)
	kappa_m = -6	(-mV)
}

STATE {
	m
}

ASSIGNED {
	ica	(mA/cm2)
	m_inf
	tadj
	eca 	(mV)
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

	m_inf = 1 / (1 + (Exp((v - theta_m)/ kappa_m))): / tadj, valido anche exp(x)

}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(Exp(x/y)-1)
	}
}

FUNCTION Exp(x) {
	if (exp(x) < 0.1) {
		Exp = 0.1 : Exp = 0
	}else{
		Exp = exp(x)
	}
} 

