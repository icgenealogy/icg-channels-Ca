:high voltage activated Ca2+ current

NEURON {
	SUFFIX Hva
	USEION ca READ cai, cao WRITE ica
	RANGE ghvabar, ica, gca, eca
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	ghvabar= 0.00034 (mho/cm2) <0,1e9>
	
	
}
STATE {
	u z
}
ASSIGNED {
	ica (mA/cm2)
	uinf zinf 
	utau (ms)
	ztau (ms)
	gca (mho/cm2)
	eca (mV)
	cai (mM)
	cao (mM)
	
}



INITIAL {
	rate(v)
	u = uinf
	z = zinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gca = ghvabar*u*u*z
	eca = 12.5 * log(cao/cai)
	ica = gca*(v-eca)
	
}

DERIVATIVE states {
	rate(v)
	u' = (uinf-u)/utau
	z' = (zinf-z)/ztau
}


UNITSOFF

PROCEDURE rate(v (mV)) {LOCAL  vu, vz, vx
	
	vx = -0.031*(v+37.1)
	vu = v+24.6
	
if (fabs(vu)<1e-04){
	   vu = vu+0.00001
	   uinf = 1/(1+exp(-(vu/11.3)))
	   utau = (1.25*(2/(exp(vx) + exp(-vx))))
}else{
	   uinf = 1/(1+exp(-(vu)/11.3))
	   utau = (1.25*(2/(exp(vx) + exp(-vx))))
}
	
	  vz = v+12.6

if (fabs(vz)<1e-04){
	  vz = vz+0.00001
	  zinf = 1/(1+exp(vz/18.9))
	  ztau = 420
}else{
	  zinf = 1/(1+exp(vz/18.9))
	  ztau = 420	
	}
}

UNITSON	










