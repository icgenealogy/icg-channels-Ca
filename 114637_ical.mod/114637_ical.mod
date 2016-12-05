TITLE Ca Channel
: High-threshold, long-lasting calcium channel
: by Ojvind Bernander 92-01-21.
: (db) added RANGE section to allow access to parameters from mech. browser
: (db) 4.1.98 modifications for CVode

NEURON {
	SUFFIX icalnew
	USEION ca READ eca WRITE ica
	GLOBAL inf
	RANGE gcabar, eca
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
:	celsius = 37	(degC)
	dt (ms)
	gcabar=.0006 (mho/cm2)
	eca = 115 (mV)
	caactvha = 25 (mV)
	caactslope = -4 (mV)
}
STATE { m }
ASSIGNED {
	ica (mA/cm2)
	inf[1]
	tau[3]
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar*m*m * (v - eca)
}

DERIVATIVE states {	: exact when v held constant
	mhn(v*1(/mV))
	m' = (inf[0] - m)/tau[0]
}

FUNCTION varss(v, i) {
	if (i==0) {
		varss = 1 / (1 + exp((v + caactvha)/(caactslope))) :Ca activation
	}
}

FUNCTION vartau(i) {
	if (i==0) {
		vartau = 2.0
	}
}	

PROCEDURE mhn(v) {
	TABLE inf, tau
	DEPEND celsius, caactvha, caactslope, dt
	FROM -100 TO 100 WITH 2000 : .1 mV steps.
	FROM i=0 TO 0 {
		tau[i] = vartau(i)
		inf[i] = varss(v,i)
	}
}


