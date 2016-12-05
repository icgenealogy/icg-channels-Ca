TITLE Calcium high-threshold current
 
COMMENT
  from Table 3 of "A branching dendritic model of a rodent CA3 pyramidal neurone." Traub RD et al. J Physiol. (1994) 
  implemented by Nikita Vladimirov <nikita.vladimirov@gmail.com>
ENDCOMMENT

NEURON {
        SUFFIX Ca
		USEION ca READ eca WRITE ica
        RANGE  gbar, g, minf
		GLOBAL Vm
} 
 
UNITS {
		(S)  = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER { 
		gbar = 1.0   		(S/cm2) 
		Vm   = -65 		(mV) : resting potential
		:Eca  = 75      (mV) : Ca reversal potential, absolute (140 mV relative to Vm)
}

ASSIGNED {
                eca (mV)
		v   (mV)
		ica  (mA/cm2)
		i   (mA/cm2)
		g   (S/cm2)
		minf
		mtau (ms) 
}

STATE { m }

BREAKPOINT {
		SOLVE states METHOD cnexp
		g = gbar * m^2
		i = g * (v - eca)
		ica = i
}

INITIAL {
		rates(v)
		m = minf
}

DERIVATIVE states {
        rates(v)
        m' = (minf - m) / mtau
}

PROCEDURE rates(v(mV)) {
		LOCAL  alpham, betam, small
        TABLE minf, mtau FROM -100 TO 50 WITH 200
		UNITSOFF
			alpham =  1.6 / ( 1 + exp(-0.072 * (v - Vm - 65)) )
			small = (v - Vm - 51.1) / 5 
			if( fabs(small) > 1e-6 ) {
				betam  =  0.02 * (v - Vm - 51.1) / ( exp( (v - Vm - 51.1) / 5 ) - 1)
			} else {
				betam  =  0.02 * 5 / ( 1 + small / 2 )
			}
			minf   = alpham / ( alpham + betam )
			mtau   = 1 / ( alpham + betam )
		UNITSON
}
