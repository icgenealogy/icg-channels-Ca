TITLE CaN.mod N-type Cav channel  
COMMENT

Original Mod File:
Original name 'nca.mod'
Santhakumar V, Aradi I, Soltesz I (2005) J Neurophysiol 93:437-53 
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=51781&file=/dentategyrusnet2005/nca.mod

Current version by A. Hanuschkin <AH, 2011> for:
Yim MY, Hanuschkin A, Wolfart J (2015) Hippocampus 25:297-308.
http://onlinelibrary.wiley.com/doi/10.1002/hipo.22373/abstract

Mod File history:
- fitted H-H parameter N-Ca from Jaffe DB, Ross WN, Lisman JE,  Lasser-Ross N, Miyakawa H, Johnston D (1994) Journal of Neurophysiology, Vol. 71 no. 3, 1065-1077
- Ca ion & L/T/N-Ca channels model of Aradi I, Holmes WR (1999) J Comput Neurosci 6:215-35
- Note that eCa is calculated during simulation by ccanl.mod. ecat, ecal values set in Santhakumar are not used in our model scripts.

ENDCOMMENT
 
UNITS {
        (mA) =		(milliamp)
        (mV) =		(millivolt)
        (uF) =		(microfarad)
	(molar) = 	(1/liter)
	(nA) = 		(nanoamp)
	(mM) = 		(millimolar)
	(um) = 		(micron)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}
 
NEURON { 
SUFFIX nca
:USEION nca READ enca WRITE inca VALENCE 2 
USEION ca READ eca WRITE ica
RANGE  gnca
RANGE gncabar
RANGE cinf, ctau, dinf, dtau, inca
GLOBAL eca
}
 
INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}
 
PARAMETER {
        v (mV) 
        celsius = 6.3 (degC)
        dt (ms) 
	gncabar = 1.0 (mho/cm2)
}
 
STATE {
	c d
}
 
ASSIGNED {
	gnca (mho/cm2)
	ica (mA/cm2)
	eca (mV)

	cinf dinf
	ctau (ms) dtau (ms) 
	cexp dexp      
} 

BREAKPOINT {
	SOLVE states
        gnca = gncabar*c*c*d
	ica = gnca*(v-eca)
}
 
UNITSOFF
 
INITIAL {
	trates(v)
	c = cinf
	d = dinf
}

PROCEDURE states() {	:Computes state variables c and d 
        trates(v)	:      at the current v and dt.
	c = c + cexp*(cinf-c)
	d = d + dexp*(dinf-d)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
LOCAL q10

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  alpha, beta, sum
        q10 = 3^((celsius - 6.3)/10) : q10=1, for 6.3 celsius
                :"c" NCa activation system
        alpha = -0.19*vtrap(v-19.88,-10)
	beta = 0.046*exp(-v/20.73)
	sum = alpha+beta        
	ctau = 1/sum      cinf = alpha/sum
                :"d" NCa inactivation system
	alpha = 0.00016*exp(-v/48.4)
	beta = 1/(exp((-v+39)/10)+1)
	sum = alpha+beta        
	dtau = 1/sum      dinf = alpha/sum
}
 
PROCEDURE trates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc
        TABLE  cinf, cexp, dinf, dexp, ctau, dtau
	DEPEND dt, celsius FROM -100 TO 100 WITH 200
                           
	rates(v)	: not consistently executed from here if usetable_hh == 1
		: so don't expect the tau values to be tracking along with
		: the inf values in hoc

        tinc = -dt * q10
	cexp = 1 - exp(tinc/ctau)
	dexp = 1 - exp(tinc/dtau)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{  
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON

