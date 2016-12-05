TITLE calcium T channel for STh

COMMENT
 Low threshold calcium channel (T-type), Wang et al. 1991 
 & Coulter et al 1989.  The original data was recorded at 22-24degC. 

 How the q10 works: There is a q10 for the rates (alpha and beta's)
 called Q10 and a Q10 for the maximum conductance called gmaxQ10.  The
 q10s should have been measured at specific temperatures temp1 and
 temp2 (that are 10degC apart). Ideally, as Q10 is temperature
 dependant, we should know these two temperatures.  We used to
 follow the more formal Arrhenius derived Q10 approach.  The
 temperature at which this channel's kinetics were recorded is tempb
 (base temperature).  What we then need to calculate is the desired
 rate scale for now working at temperature celsius (rate_k).  This was
 given by the empirical Arrhenius equation, using the Q10, but now is 
 using the quick Q10 approximation. 
ENDCOMMENT

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
     FARADAY = (faraday) (coulomb)
           R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX CaT
	USEION ca READ cai,cao,eca WRITE ica
	RANGE gcaT, iCaT
	GLOBAL activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gcaT  = 0.001 (mho/cm2)
	iCaT  = 0.0 (mA/cm2)
	eca
	cai
	cao
	celsius

	activate_Q10 = 1
	Q10 = 1.515804730e+00
	gmaxQ10 = 1.515804730e+00
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

STATE {
        r s d
}

ASSIGNED { 
        ica (mA/cm2)
	ralpha (/ms)
	rbeta (/ms)
	salpha (/ms)
	sbeta (/ms)
	dalpha (/ms)
	dbeta (/ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica  = (gcaT*gmax_k)*r*r*r*s*ghkg(v,cai,cao,2)
	iCaT = ica
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
	  rate_k = Q10^((celsius-tempb)/10)
          gmax_k = gmaxQ10^((celsius-tempb)/10)
	}else{
	  rate_k = 1.0
	  gmax_k = 1.0
	}
        settables(v)
        r = ralpha/(ralpha+rbeta)
        s = (salpha*(dbeta+dalpha) - (salpha*dbeta))/((salpha+sbeta)*(dalpha+dbeta) - (salpha*dbeta))
	d = (dbeta*(salpha+sbeta) - (salpha*dbeta))/((salpha+sbeta)*(dalpha+dbeta) - (salpha*dbeta))
}

DERIVATIVE states {  
	settables(v)      :Computes state variables at the current v and dt.
	r' = ((ralpha*(1-r)) - (rbeta*r))
	d' = ((dbeta*(1-s-d)) - (dalpha*d))
	s' = ((salpha*(1-s-d)) - (sbeta*s))
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
			  :Voltage shift (for temp effects) of -1.9278 added
        LOCAL   bd
        TABLE ralpha, rbeta, salpha, sbeta, dalpha, dbeta DEPEND celsius FROM -100 TO 100 WITH 400

		:"r" CaT activation system
	ralpha = rate_k * 1.0/(1.7+exp(-(v + 26.2722)/13.5))
	rbeta  = rate_k * exp(-(v + 61.0722)/7.8)/(exp(-(v + 26.8722)/13.1)+1.7)

                :"s" CaT fast inactivation system
        salpha = rate_k * exp(-(v + 158.3722)/17.8)
        sbeta  = rate_k * (sqrt(0.25+exp((v + 81.5722)/6.3))-0.5) * (exp(-(v + 158.3722)/17.8))

	        :"d" CaT slow inactivation system
	bd     = sqrt(0.25+exp((v + 81.5722)/6.3))
	dalpha = rate_k * (1.0+exp((v + 35.4722)/30.0))/(240.0*(0.5+bd))
        dbeta  = rate_k * (bd-0.5)*dalpha
}

UNITSON

INCLUDE "ghk.inc"
