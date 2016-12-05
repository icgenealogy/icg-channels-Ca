NEURON {
	SUFFIX cap_gp
	USEION ca READ cai, cao WRITE ica
	RANGE gbar, ica
	GLOBAL minf,mtau
	GLOBAL monovalConc, monovalPerm
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(mM) = (milli/liter)
	(S) = (siemens)
	F = 9.6485e4   (coul)
	R = 8.3145 (joule/degC)
	(mC) = (millicoul)
}

PARAMETER {
	v			(mV)

	gbar = 1		(cm/s)
	monovalConc = 140	(mM)
	monovalPerm = 0

	cai			(mM)
	cao			(mM)
	celsius (degC)
}

ASSIGNED {
	ica	(mA/cm2)
        minf
	mtau	(ms)
	T	(degC)
	E	(volts)
	g 	(cm/s)
}

STATE {
	m
}

INITIAL {
	rates(v)
	m = minf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m
	ica = g * ghk(v, cai, cao, 2)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/mtau
}

FUNCTION ghk( v(mV), ci(mM), co(mM), z)  (mC/cm3) { LOCAL Ci
	T = 22 + 273.19  : Kelvin
        E = v * 1e-3 (volts/mV)
        Ci = ci + (monovalPerm) * (monovalConc)        : Monovalent permeability
	if (fabs(1-exp(-z*(F*E)/(R*T))) < 1e-6) { : denominator is small -> Taylor series
		ghk = (1e-3) * z * F * (Ci-co*exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
	} else {
		ghk = (1e-3) * z^2*(E*F^2)/(R*T)*(Ci-co*exp(-z*(F*E)/(R*T)))/(1-exp(-z*(F*E)/(R*T)))
	}
}

PROCEDURE rates (v (mV)) {
	LOCAL q10
	minf = 1/(1+exp(-(v + 19 (mV) ) / 5.5 (mV)))
	q10 = 3^((celsius - 22 (degC))/10 (degC) )
	mtau = q10*(mtau_func(v)) * 1e3
}

FUNCTION mtau_func( v (mV) ) (ms) {
        if (v > -50) {
            mtau_func = .000191 (ms) + .00376 (ms) *exp(-((v+41.9 (mV))/27.8 (mV) )^2)
        } else {
            mtau_func = .00026367 (ms) + .1278 (ms) * exp(v* .10327 (1/mV))
        }
}
