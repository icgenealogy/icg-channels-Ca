COMMENT
	high-threshold calcium channel from Av-Ron and Vidal, 1999
	Implemented by C. Weaver, 2003
ENDCOMMENT

UNITS {
	(molar) =	(1/liter)
	(mM) =	(millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)

}

NEURON {
	SUFFIX cahi
	USEION ca READ eca WRITE ica
        RANGE gbar
        GLOBAL xinf
	RANGE tot
}

PARAMETER {
	v (mV)
	celsius 	(degC)
	: gbar=.001 (mho/cm2)
	gbar = 0.0002 (mho/cm2)
	: cai (mM)
	: cao (mM)
	xtau=5 (ms)
	Kc=1   (mM)
	ax=0.08 (/mV)
	vhx=-30 (mV)
	vrest = 124	(mV)
	simp = 0
}


STATE {
	x
}

ASSIGNED {
	ica (mA/cm2)
	tot (mA/cm2)
	cai (mM)
        gca (mho/cm2)
        xinf
	eca (mV)
}

INITIAL {
	rate(v)
	x = xinf
: printf("cahi ica=%g\n",ica)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gca = gbar*x*x*Kc/(Kc+cai)
	tot = gca*(v-eca)
	if( simp > 0 ) {
		: Av-Ron 1991, simpler
		gca = gbar*x
		: gca = gbar*x*x
	}
	: ica = gca*(v-vrest)
	ica = gca*(v-eca)
}

FUNCTION expn(v (mV),a(/mV), vhalf(mV)) {
  	expn = exp(-2*a*(v-vhalf))
}

DERIVATIVE state {  
        rate(v)
        x' = (xinf - x)/xtau
}

PROCEDURE rate(v (mV)) { :callable from hoc
	xinf = 1/(1 + expn(v,ax,vhx))
}
 







