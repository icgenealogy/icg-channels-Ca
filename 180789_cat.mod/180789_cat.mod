TITLE T calcium current
: FORREST MD (2014) Two Compartment Model of the Cerebellar Purkinje Neuron
 
COMMENT
  from "An Active Membrane Model of the Cerebellar Purkinje Cell
        1. Simulation of Current Clamp in Slice"
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX cat
        USEION ca READ cai, cao WRITE ica
        RANGE  gcabar, ica, gca, minf, hinf, mexp, hexp
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gcabar = .0005 (mho/cm2)
        eca = 135 (mV)
	cai	= 0.40e-4 (mM)		: adjusted for eca=135mV
	cao	= 2.4	(mM)

}
 
STATE {
        m h
}
 
ASSIGNED {
        ica (mA/cm2)
        gca minf hinf mexp hexp 
}
 
BREAKPOINT {
        SOLVE states
        gca = gcabar * m*h
	ica = gca* (v-eca)
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {  :Computes state variables m, h
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                :"m" calcium activation system
        alpha = 2.6/(1+exp((v+21)/(-8)))
        beta =  0.18/(1+exp((v+40)/4))
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
                :"h" calcium inactivation system
        alpha = 0.0025/(1+exp((v+40)/8))
        beta = 0.19/(1+exp((v+50)/(-10)))
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
}

 
UNITSON

