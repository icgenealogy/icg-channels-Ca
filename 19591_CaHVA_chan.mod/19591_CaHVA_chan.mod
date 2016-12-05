TITLE CaHVA_chan.mod  High voltage activated ca channel, adapted from Traub (1991). 
 
COMMENT
%W%                                 %G%
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX cahva_chan
        USEION ca READ eca WRITE ica
        RANGE gbar, i
        GLOBAL minf, hinf,
               a1m, b1m, c1m, a2m, b2m, c2m, taum_min,
               a1h, b1h, c1h, d1h, e1h, a2h, b2h, c2h, d2h, e2h, tauh_min
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        gbar = 2.9089e-3 (mho/cm2)
        eca (mV)

        a1m = 1.6
        b1m = -0.072
        c1m = -5
        a2m = 0.02
        b2m = 0.2
        c2m = 8.9
        taum_min = 1e-3

        a1h = 0.005
        b1h = 0
        c1h = 0.005
        d1h = -0.05
        e1h = -60
        a2h = 0
        b2h = 0.005
        c2h = -0.005
        d2h = -0.05
        e2h = -60
        tauh_min = 1e-3
}
 
STATE {
        m h
}
 
ASSIGNED {
        ica (mA/cm2)
        minf hinf
        i
}
 
LOCAL mexp, hexp
 
BREAKPOINT {
        SOLVE states
        i  = gbar*m*m*h*(v - eca)
        ica = i
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
}

PROCEDURE states() {  :Computes state variables m, h, and n 
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        VERBATIM
        return 0;
        ENDVERBATIM
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  tau,alpha,beta
        TABLE minf, mexp, hinf, hexp DEPEND dt FROM -100 TO 100 WITH 2000

                :"m" HVA activation system
        alpha = a1m/(1+exp(b1m*(v+c1m)))
        beta = a2m*(v+c2m)/(exp(b2m*(v+c2m))-1)
        tau = 1/(alpha + beta)
        minf = alpha*tau
        if (tau<taum_min) { tau = taum_min }
        mexp = 1 - exp(-dt/tau)

                :"h" HVA inactivation system
        if (v<e1h) {
          alpha = a1h
        }
        else {
          alpha = b1h+c1h*exp(d1h*(v-e1h))
        }
        if (v<e2h) {
          beta = a2h
        }
        else {
          beta = b2h+c2h*exp(d2h*(v-e2h))
        }
        tau = 1/(alpha + beta)
        hinf = alpha*tau
        if (tau<tauh_min) { tau = tauh_min }
        hexp = 1 - exp(-dt/tau)
}
 
 
UNITSON

