TITLE CaT channel alpha-1H from McRory et al, 2001
: Reversal potential described by Nernst equation
: M.Migliore Jan 2003

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v (mV)
	celsius		(degC)
	gbar=.008 (mho/cm2)
        vhalfn=-43.15   (mV)
        vhalfl=-73.9   (mV)
        kn=5.34   (1)
        kl=-2.76   (1)
	q10=2.3
	cai 	= .00005 (mM)	: initial [Ca]i = 50 nM
	cao 	= 2	(mM)	: [Ca]o = 2 mM
	eca
}


NEURON {
	SUFFIX cat1h
	USEION ca READ cai,cao,eca WRITE ica
        RANGE gbar, carev
        GLOBAL ninf,linf,taul,taun, q10
}

STATE {
	n
        l
}

ASSIGNED {
	ica	(mA/cm2)		: current
	carev	(mV)			: rev potential
        ninf
        linf      
        taul
        taun
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gbar*n*l*(v-carev)
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-22)/10)
        ninf = 1/(1 + exp(-(v-vhalfn)/kn))
        linf = 1/(1 + exp(-(v-vhalfl)/kl))
        taun = (0.774+0.14*exp(-v/13.27))/qt
        taul = (22.25+0.0455*exp(-v/7.46))/qt
}














