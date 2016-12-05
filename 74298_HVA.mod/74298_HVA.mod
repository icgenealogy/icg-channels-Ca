TITLE calcium HVA channels for STh

COMMENT
 High threshold calcium channel (N/L-type), Brown et al. 1993. & Fox
 et al. (1989).  Both done at temperature 22degC.

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

 Adding CaL [Ca]i dependent inactivation.  This is only for the L-type
 component, and is called inactivation variable 'h'.  
ENDCOMMENT

UNITS {
	(mM) = (milli/liter)
	(mV) = (millivolt)
	(mA) = (milliamp)
     FARADAY = (faraday) (coulomb)
           R = (k-mole) (joule/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX HVA
	USEION ca READ cai,cao,eca WRITE ica
	RANGE gcaN, gcaL, iNCa, iLCa
	GLOBAL inactLtau,inactLmax,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
        gcaL  = 0.002 (mho/cm2)
	gcaN  = 0.012 (mho/cm2)
	iNCa  = 0.0 (mA/cm2)
	iLCa  = 0.0 (mA/cm2)
	inactLtau = 1220.0 (ms)
	inactLmax = 5.291291201e-01
	eca
	cai
	cao
	celsius

	activate_Q10 = 1
	Q10 = 1.948259241e+00
	gmaxQ10 = 1.948259241e+00
	temp1 = 20.0 (degC)
	temp2 = 30.0 (degC)
	tempb = 22.0 (degC)
}

STATE {
        q u h
}

ASSIGNED { 
        ica (mA/cm2)
	qinf
	uinf
	hinf
	qtau (ms)
	utau (ms)
	htau (ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	LOCAL vghk
	SOLVE states METHOD cnexp
	vghk = ghkg(v,cai,cao,2)
	iNCa = gmax_k*(gcaN * u)*q*q*vghk
	iLCa = gmax_k*(gcaL)*q*q*h*vghk
	ica  = iNCa + iLCa
}


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
	q = qinf
	u = uinf
	setCadepLinact(cai)
	h = hinf
}

DERIVATIVE states {  
	settables(v)  
	q' = (qinf-q)/qtau
	u' = (uinf-u)/utau
	setCadepLinact(cai)
	h' = (hinf-h)/htau
}

PROCEDURE settables(v) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
			  :Voltage shifts (for temp effects) of -8.25 and -14.67 added respt.
        TABLE qinf, qtau, uinf, utau DEPEND celsius FROM -100 TO 100 WITH 400

                :"q" N/L Ca activation system
        qinf   = 1.0/(1.0 + exp((-16.3547869 - v)/11.3))
        qtau   = (1.25/(cosh(-0.031 * (v + 28.8547869)))) /rate_k

                :"u" N inactivation system - voltage dependent.
        uinf   = 1.0/(1.0 + exp((v + 45.3326653)/12.5))
        utau   = (98.0 + cosh(0.021*(24.7673347-v))) /rate_k
}

PROCEDURE setCadepLinact(cai) { : set Ca dependent L-type calcium channel inactivation
                :"h" L inactivation system - [Ca]i dependent.
	hinf   = inactLmax+((1.0-inactLmax)/(1.0 + exp((cai-0.7)/0.15)))
        htau   = inactLtau /rate_k
}

INCLUDE "ghk.inc"


