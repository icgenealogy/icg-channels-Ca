
COMMENT

ca.mod
Uses fixed eca instead of GHK eqn

HVA Ca current
Based on Reuveni, Friedman, Amitai and Gutnick (1993) J. Neurosci. 13:
4609-4621.

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ca
	USEION ca READ eca WRITE ica
	RANGE m, h, gca, gbar
	RANGE minf, hinf, mtau, htau
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
	gbar = 0.1   	(pS/um2)	: 0.12 mho/cm2
	vshift = 0	(mV)		: voltage shift (affects all)

	cao  = 2.5	(mM)	        : external ca concentration
	cai		(mM)
						
	temp = 23	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	PI	= (pi) (1)
	(molar) = (1/liter)
	(mM) = (millimolar)
} 

ASSIGNED {
	ica 		(mA/cm2)
	gca		(pS/um2)
	eca		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
:printf("ca INTIIAL entry t=%g cai=%g ica=%g eca=%g\n", t, cai, ica, eca)
	trates(v+vshift)
	m = minf
	h = hinf
:printf("ca INTIIAL exit t=%g cai=%g ica=%g eca=%g\n", t, cai, ica, eca)
}

BREAKPOINT {
        SOLVE states METHOD cnexp
:printf("ca BREAKPOINT entry t=%g cai=%g ica=%g eca=%g\n", t, cai, ica, eca)
        gca = tadj*gbar*m*m*h
	ica = (1e-4) * gca * (v - eca)
:printf("ca BREAKPOINT exit t=%g cai=%g ica=%g eca=%g\n", t, cai, ica, eca)
} 


DERIVATIVE states {
:printf("ca states entry t=%g cai=%g ica=%g eca=%g\n", t, cai, ica, eca)
        trates(v+vshift)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
:printf("ca states exit t=%g cai=%g ica=%g eca=%g\n", t, cai, ica, eca)
}


PROCEDURE trates(v (mV)) {  
                      
        TABLE minf, mtau, hinf, htau
	DEPEND celsius, temp
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

UNITSOFF
        tadj = q10^((celsius - temp)/10)
UNITSON
        mtau = mtau/tadj
        htau = htau/tadj
}


UNITSOFF
PROCEDURE rates(vm (mV)) {
        LOCAL  a, b

	a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1)
	b = 0.94*exp((-75-vm)/17)
	
	mtau = 1/(a+b)
	minf = a*mtau

		:"h" inactivation 

	a = 0.000457*exp((-13-vm)/50)
	b = 0.0065/(exp((-vm-15)/28) + 1)

	htau = 1/(a+b)
	hinf = a*htau
}
UNITSON

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
