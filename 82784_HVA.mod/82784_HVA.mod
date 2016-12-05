: High voltage activated Ca2+ current
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
	SUFFIX HVA
	USEION ca READ cai, cao WRITE ica
	RANGE gHVAbar, gca, eca, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gHVAbar= 0.00034 (mho/cm2) <0,1e9>
}

ASSIGNED {
	v    (mV)
	ica (mA/cm2)
	uinf zinf 
	utau (ms)
	ztau (ms)
	gca (mho/cm2)
	eca (mV)
	cai (mM)
	cao (mM)
}

STATE {
	u 
        z 
}

INITIAL {
	rate(v)
	u = uinf
	z = zinf
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	gca = gHVAbar*u*u*z
	eca = 12.5*log(cao/cai)
	ica = gca*(v-eca)
}

DERIVATIVE states {
	rate(v)
	u' = (uinf-u)/utau
	z' = (zinf-z)/ztau
}

UNITSOFF

PROCEDURE rate(v (mV)) {
	LOCAL vx
	vx = -0.031*(v+37.1)
	uinf = 1/(1+exp(-(v+24.6)/11.3))
	utau = 1.25*(2/(exp(vx)+exp(-vx)))
	zinf = 1/(1+exp((v+12.6)/18.9))
	ztau = 420/3
}

UNITSON	
