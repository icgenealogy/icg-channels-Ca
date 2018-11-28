
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
	SUFFIX casoma
	USEION ca READ eca WRITE ica
	RANGE m, h, gca, gbar
	RANGE minf, hinf, mtau, htau
	GLOBAL q10, temp, tadj, vmin, vmax, vshift
	
	RANGE v1, v2, v3, v4, v5, v6, v7, v8
}

PARAMETER {
	gbar   = 1.25(pS/um2)	: 0.12 mho/cm2
	vshift = 0	(mV)		: voltage shift (affects all)

	cao  = 2.5	(mM)	        : external ca concentration
	cai		(mM)
						
	temp = 30	(degC)		: original temp (23 Jary Y. Delgado)
	q10  = 2.3			  : temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	
	v1 = 13 (mV)
	v2 = 50 (mV)
	v3 = 15 (mV)
	v4 = 8.5 (mV)
	
	v5 = 27 (mV)
	v6 = 6 (mV)
	v7 = 90 (mV)
	v8 = 8.5 (mV)



	
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
      (molar) = (1/liter)
      (mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
	PI	= (pi) (1)
} 

ASSIGNED {
	ica 		(mA/cm2)
	gca		(pS/um2)
	eca		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
	
	a
	b
}
 

STATE { m h }

INITIAL { 
	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states
        gca = tadj*gbar*m*m*h
	  ica = (1e-4) * gca * (v - eca)
} 

LOCAL mexp, hexp

PROCEDURE states() {
     
	 
	 :   trates(v+vshift)  

	 LOCAL tinc
	 
:	a = 0.055*(-27 - v)/(exp((-27-v)/3.8) - 1)
:	b = 0.94*exp((-75-v)/17)

	a = 0.055*(-v5 - v)/(exp((-v5-v)/v6) - 1)
	b = 0.94*exp((-v7-v)/v8)

	
	mtau = 1/(a+b)
	minf = a*mtau

	:"h" inactivation

:	a = 0.000457*exp((-13-v)/50)
:	b = 0.0065/(exp((-v-15)/28) + 1)

	a = 0.000457*exp((-v1-v)/v2)
	b = 0.0065/(exp((-v-v3)/v4) + 1)

	
	htau = 1/(a+b)
	hinf = a*htau

	tadj = q10^((celsius - temp)/10)
    tinc = -dt * tadj

    mexp = 1 - exp(tinc/mtau)
    hexp = 1 - exp(tinc/htau)
	
	
		
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
	VERBATIM
	return 0;
	ENDVERBATIM
}


PROCEDURE trates(v (mV)) {  
                      
        LOCAL tinc
        TABLE minf, mexp, hinf, hexp
	DEPEND dt, celsius, temp
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

        tadj = q10^((celsius - temp)/10(degC))
        tinc = -dt * tadj

        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
}


PROCEDURE rates(vm (mV)) {  
        LOCAL  a, b

	a = 0.055(/ms/mV)*(-27(mV) - vm)/(exp((-27(mV)-vm)/3.8(mV)) - 1)
	b = 0.94(/ms)*exp((-75(mV)-vm)/17(mV))
	
	mtau = 1/(a+b)
	minf = a*mtau

	:"h" inactivation

	a = 0.000457(/ms)*exp((-13(mV)-vm)/50(mV))
	b = 0.0065(/ms)/(exp((-vm-15(mV))/28(mV)) + 1)

	
	htau = 1/(a+b)
	hinf = a*htau
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
