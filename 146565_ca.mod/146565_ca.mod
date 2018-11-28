COMMENT
https://senselab.med.yale.edu/modeldb/ShowModel.asp?model=2488&file=\patdemo\ca.mod

ca.mod
Uses fixed eca instead of GHK eqn

HVA Ca current
Based on Reuveni, Friedman, Amitai and Gutnick (1993)
J. Neurosci. 13:4609-4621.

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu

26 Ago 2002 Modification of original channel to allow 
variable time step and to correct an initialization error.
Done by Michael Hines(michael.hines@yale.e) and 
Ruggero Scorcioni(rscorcio@gmu.edu) at EU Advance Course 
in Computational Neuroscience. Obidos, Portugal

20110202 made threadsafe by Ted Carnevale

Special comment:

This mechanism was designed to be run at a single operating 
temperature--37 deg C--which can be specified by the hoc 
assignment statement
celsius = 37
This mechanism is not intended to be used at other temperatures, 
or to investigate the effects of temperature changes.

Zach Mainen created this particular model by adapting conductances 
from lower temperature to run at higher temperature, and found it 
necessary to reduce the temperature sensitivity of spike amplitude 
and time course.  He accomplished this by increasing the net ionic 
conductance through the heuristic of changing the standard HH 
formula
  g = gbar*product_of_gating_variables
to
  g = tadj*gbar*product_of_gating_variables
where
  tadj = q10^((celsius - temp)/10)
  temp is the "reference temperature" (at which the gating variable
    time constants were originally determined)
  celsius is the "operating temperature"

Users should note that this is equivalent to changing the channel 
density from gbar at the "reference temperature" temp (the 
temperature at which the at which the gating variable time 
constants were originally determined) to tadj*gbar at the 
"operating temperature" celsius.
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    THREADSAFE
	SUFFIX cah
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
:	dt		(ms)
:	celsius		(degC)
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
    tadj = q10^((37 - temp)/10) : make all threads calculate tadj at initialization

	trates(v+vshift)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = tadj*gbar*m*m*h
	ica = (1e-4) * gca * (v - eca)
} 

LOCAL mexp, hexp

:PROCEDURE states() {
:        trates(v+vshift)      
:        m = m + mexp*(minf-m)
:        h = h + hexp*(hinf-h)
:	VERBATIM
:	return 0;
:	ENDVERBATIM
:}

DERIVATIVE states {
        trates(v+vshift)      
        m' =  (minf-m)/mtau
        h' =  (hinf-h)/htau
}

PROCEDURE trates(v) {  
                      
        
        TABLE minf, hinf, mtau, htau 
	DEPEND  celsius, temp
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

:        tinc = -dt * tadj

:        mexp = 1 - exp(tinc/mtau)
:        hexp = 1 - exp(tinc/htau)
}


PROCEDURE rates(vm) {  
        LOCAL  a, b

        tadj = q10^((37 - temp)/10)

	a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1)
	b = 0.94*exp((-75-vm)/17)

	mtau = 1/tadj/(a+b)
	minf = a/(a+b)

		:"h" inactivation 

	a = 0.000457*exp((-13-vm)/50)
	b = 0.0065/(exp((-vm-15)/28) + 1)

	htau = 1/tadj/(a+b)
	hinf = a/(a+b)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}