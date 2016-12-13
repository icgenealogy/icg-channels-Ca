TITLE CaL.mod L-type Cav channel
COMMENT

Mod file by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Mod File history:
- fitted H-H parameter N-Ca from Jaffe DB, Ross WN, Lisman JE,  Lasser-Ross N, Miyakawa H, Johnston D (1994) Journal of Neurophysiology, Vol. 71 no. 3, 1065-1077
- Ca ion & L/T/N-Ca channels model of  Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35
- checked and adapted by Hanuschkin in 2011
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
SUFFIX lca
:USEION lca READ elca WRITE ilca VALENCE 2 
USEION ca READ eca WRITE ica
RANGE  glca
RANGE glcabar
RANGE einf, etau, ica
GLOBAL eca
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
	glcabar = 1.0 (mho/cm2)
}
 
STATE {
	e
}
 
ASSIGNED {
	glca (mho/cm2)
	ica (mA/cm2)
	eca (mV)

	einf 
	etau (ms) 
	eexp      
} 

BREAKPOINT {
	SOLVE states
        glca = glcabar*e*e
	ica = glca*(v-eca)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	e = einf
}

PROCEDURE states() {	:Computes state variables e 
        trates(v)	:      at the current v and dt.
	e = e + eexp*(einf-e)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10) : q10=1 for 6.3 celcius
                :"e" LCa activation system
        alpha = -15.69*vtrap(v-81.5,-10)	
	beta = 0.29*exp(-v/10.86)	
	sum = alpha+beta        
	etau = 1/sum      einf = alpha/sum
                :no LCa inactivation system
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
        TABLE  einf, eexp, etau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

	tinc = -dt * q10
	eexp = 1 - exp(tinc/etau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

