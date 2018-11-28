: Ca channels (T,N,L-type)
: Aradi and Holmes (1999)

NEURON {
	SUFFIX Ca
	USEION ca READ eca WRITE ica
	RANGE gtcabar, gncabar, glcabar, gtca, gnca, glca
	GLOBAL ca0, cao
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	B = .26 (mM-cm2/mA-ms)
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	TEMP = 25 (degC)
}

PARAMETER {
	ca0 = .00007	(mM)		: initial calcium concentration inside
	cao = 2		(mM)		: calcium concentration outside
	tau = 9		(ms)
	gtcabar = .01	(S/cm2)	: maximum permeability
	gncabar = .01	(S/cm2)
	glcabar = .01	(S/cm2)
}

ASSIGNED {
	v		(mV)
	eca		(mV)
	ica		(mA/cm2)
	gtca		(S/cm2)
	gnca		(S/cm2)
	glca		(S/cm2)
	i		(nA)
}

STATE { ca_i (mM) a b c d e}

BREAKPOINT {
	SOLVE state METHOD cnexp
	:e_ca = (1000)*(TEMP+273.15)*R/(2*F)*log(cao/ca_i)
	gtca = gtcabar*a*a*b
	gnca = gncabar*c*c*d
	glca = glcabar*e*e
	ica = (gtca+gnca+glca)*(v - eca)   : only L type

}

DERIVATIVE state {	: exact when v held constant; integrates over dt step
	ca_i' = -B*ica-(ca_i-ca0)/tau
	a' = alphaa(v)*(1-a)-betaa(v)*a
	b' = alphab(v)*(1-b)-betab(v)*b
	c' = alphac(v)*(1-c)-betac(v)*c
	d' = alphad(v)*(1-d)-betad(v)*d
	e' = alphae(v)*(1-e)-betae(v)*e
}

INITIAL {
	ca_i = ca0
	a = alphaa(v)/(alphaa(v)+betaa(v))
	b = alphab(v)/(alphab(v)+betab(v))
	c = alphac(v)/(alphac(v)+betac(v))
	d = alphad(v)/(alphad(v)+betad(v))
	e = alphae(v)/(alphae(v)+betae(v))
}

FUNCTION alphaa(v (mV)) (/ms) {
	alphaa = f(2,0.1,v,19.26)
}

FUNCTION betaa(v (mV)) (/ms) {
	betaa = exponential(0.009,-0.045393,v,0)
}

FUNCTION alphab(v (mV)) (/ms) {
	alphab = exponential(1e-6,-0.061501,v,0)
}

FUNCTION betab(v (mV)) (/ms) {
	betab = logistic(1,-0.1,v,29.79)
}

FUNCTION alphac(v (mV)) (/ms) {
	alphac = f(1.9,0.1,v,19.88)
}

FUNCTION betac(v (mV)) (/ms) {
	betac = exponential(0.046,-0.048239,v,0)
}

FUNCTION alphad(v (mV)) (/ms) {
	alphad = exponential(1.6e-4,-0.020661,v,0)
}

FUNCTION betad(v (mV)) (/ms) {
	betad = logistic(1,-0.1,v,39)
}

FUNCTION alphae(v (mV)) (/ms) {
	alphae = f(156.9,0.1,v,81.5)
}

FUNCTION betae(v (mV)) (/ms) {
	betae = exponential(0.29,-0.092081,v,0)
}

FUNCTION f(A, k, v (mV), D) (/ms) {
	LOCAL x
	UNITSOFF
	x = k*(v-D)
	if (fabs(x) > 1e-6) {
		f = A*x/(1-exp(-x))
	}else{
		f = A/(1-0.5*x)
	}
	UNITSON
}

FUNCTION logistic(A, k, v (mV), D) (/ms) {
	UNITSOFF
	logistic = A/(1+exp(k*(v-D)))
	UNITSON
}

FUNCTION exponential(A, k, v (mV), D) (/ms) {
	UNITSOFF
	exponential = A*exp(k*(v-D))
	UNITSON
}
