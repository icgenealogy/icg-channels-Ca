TITLE ICa.mod    

COMMENT
A very simple Ca2+ channel

Author: Fredrik Edin, 2003
Address: freedin@nada.kth.se

ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

NEURON {
        SUFFIX ICa
        USEION ca READ eca WRITE ica
        RANGE gcabar, gca, m
}
 
PARAMETER {
        gcabar 	= 0.0015 (mho/cm2)	<0,1e9>	
	:eca	= 120 	 (mV)
}
 
ASSIGNED {
        eca (mV)
        v 	(mV)
	gca 	(mho/cm2)
        ica 	(mA/cm2)
	m	(1)
}
 
BREAKPOINT {
	UNITSOFF
        m =  1 / ( 1 +  exp(-(v+20)/9) )
	UNITSON
        gca = gcabar*m*m
	ica = gca * ( v - eca )
}

