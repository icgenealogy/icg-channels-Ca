TITLE CA1CaG.mod  - generic HVA Ca
 
COMMENT
From Warman, Durand and Yuen J. Neurophys. 71:2033-2045, 1994
Based on Kay and Wong (1987) data.
As used by Lipowsky et al (1996) with fixed ECa
BPG 29-10-99
Scaling of inactivation time constant (tc) added
BPG 5-1-01
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CA1CaG
        USEION ca READ eca WRITE ica
        RANGE gcabar,gca
        GLOBAL minf, hinf, mexp, hexp, tc
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 36 (degC)
        dt (ms)
        gcabar = 0.01 (mho/cm2)
        :eca = 80 (mV)
        tc = 1 (1)
}
 
STATE {
        m h
}
 
ASSIGNED {
        eca (mV)
        ica (mA/cm2)
        minf hinf mexp hexp
}
 
BREAKPOINT {
        SOLVE states
        ica = gcabar*m*m*h*(v - eca)
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
        LOCAL  alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, tc FROM -100 TO 100 WITH 200
    :"m" calcium activation system
        alpha = -0.16 * vtrap(v+26,-4.5)
        beta =  0.04 * vtrap(v+12,10)
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(-dt*sum)
    :"h" calcium inactivation system
        alpha = 2 / exp((v+94)/10)
        beta = 8 / (exp(-(v-68)/27) + 1)
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(-dt*sum/tc)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
