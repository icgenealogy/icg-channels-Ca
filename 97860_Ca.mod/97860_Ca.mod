: HH-style Ca channel using parameters from Table 1 in
: Wilson and Callaway, 2000, J Neurophysiology

NEURON {
	SUFFIX Ca
	USEION ca READ eca,cai WRITE ica
	RANGE g, gbar
	GLOBAL vh, vc
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER {
	gbar = 1	(S/cm2)
	vh = -40	(mV)
	vc = 7		(mV)
	:eca = 100	(mV)
}

ASSIGNED {
        eca (mV)
	g       (S/cm2)
	v		(mV)
	ica		(mA/cm2)
	cai		(mM)
}

BREAKPOINT {
	values()
	ica = g*(v-eca)
}

INITIAL {
	values()
}

PROCEDURE values() {
	g = gbar/(1 + exp(-(v-vh)/vc))
}
