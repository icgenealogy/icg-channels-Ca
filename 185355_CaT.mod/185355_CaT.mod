TITLE CaT.mod T-type Cav channel
COMMENT

Mod File by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Mod File history:
- fitted H-H parameter N-Ca from Jaffe DB, Ross WN, Lisman JE,  Lasser-Ross N, Miyakawa H, Johnston D (1994) Journal of Neurophysiology, Vol. 71 no. 3, 1065-1077
- Ca ion & L/T/N-Ca channels model of  Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35
- Note that eCa is calculated during simulation by ccanl.mod. ecat, ecal values set in Santhakumar are not used in our model scripts.

ENDCOMMENT
 
UNITS {
        (mA) =		(milliamp)
        (mV) =		(millivolt)
        (uF) = 		(microfarad)
	(molar) = 	(1/liter)
	(nA) = 		(nanoamp)
	(mM) = 		(millimolar)
	(um) = 		(micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
SUFFIX tca
:USEION tca READ etca WRITE itca VALENCE 2 
USEION ca READ eca WRITE ica
RANGE gtca
RANGE gcatbar
RANGE ainf, atau, binf, btau, itca
GLOBAL eca
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
	gcatbar = 1.0 (mho/cm2)
}
 
STATE {
	a b
}
 
ASSIGNED {
        gtca (mho/cm2)
	ica (mA/cm2)
	eca (mV)

	ainf binf
	atau (ms) btau (ms) 
	aexp bexp      
} 

BREAKPOINT {
	SOLVE states
        gtca = gcatbar*a*a*b
	ica = gtca*(v-eca)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	a = ainf
	b = binf
}

PROCEDURE states() {	:Computes state variables a and b 
        trates(v)	:      at the current v and dt.
	a = a + aexp*(ainf-a) : i.e. a_{t+1} = a_t*exp(-dt/atau)+ainf*(1-exp(-dt/atau)); da/dt = 1/atau*(ainf-a)
	b = b + bexp*(binf-b)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10) : q10 = 1 for 6.3 celsius
                :"a" TCa activation system
        alpha = -0.2*vtrap(v-19.26,-10)		
	beta = 0.009*exp(-v/22.03)		
	sum = alpha+beta        
	atau = 1/sum      ainf = alpha/sum
                :"b" TCa inactivation system
	alpha = 1e-6*exp(-v/16.26)		
	beta = 1/(exp((29.79-v)/10)+1)		
	sum = alpha+beta        
	btau = 1/sum      binf = alpha/sum
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
        TABLE  ainf, aexp, binf, bexp, atau, btau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

	       tinc = -dt * q10
	aexp = 1 - exp(tinc/atau)
	bexp = 1 - exp(tinc/btau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

