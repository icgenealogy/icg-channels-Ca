TITLE calcium HVA channels for GPi neuron model

COMMENT

 High threshold calcium channel (N/L-type), Brown et al. 1993. & Fox
 et al. (1989).  Both done at temperature 22degC.  Implemented from
 Gillies2006.

 Adding CaL [Ca]i dependent inactivation.  This is only for the L-type
 component, and is called inactivation variable 'h'.  

 Q10=1.95 --> rate_k=exp(log(Q10)*((1/295)-(1/309))/((1/293)-(1/303)))=2.49

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
    GLOBAL inactLtau,inactLmax,rate_k,gmax_k
}

PARAMETER {
    v (mV)
    dt (ms)
    gcaL  = 0.002 (mho/cm2)
    gcaN  = 0.012 (mho/cm2)
    iNCa  = 0.0 (mA/cm2)
    iLCa  = 0.0 (mA/cm2)
    inactLtau = 1220.0 (ms)
    inactLmax = 0.529
    eca
    cai
    cao
    celsius
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
    rate_k = 2.49
    gmax_k = 2.49
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

:INCLUDE "custom_code/inc_files/114685_ghk.inc"
: included here for auto-launch on different platforms until 
: finding file in path issue resolved
FUNCTION ghkg(v(mV), ci(mM), co(mM), z) (mV) {
    
    : Called by Ih.mod

    LOCAL nu,ff,enu,fnu
    : Here we calculate an effective drive from the GHK equation
    : define
    :    f   = 10^3 RT/(zF)
    :    nu  = v/f  
    :        = z v10^-3 F / (RT) 
    : note the 10e-3 converts [mV] to [V]
    :    nu  = z V F / (RT)
    :
    :    enu = exp(nu)
    :        = exp(z V F / (RT))
    :
    :    fnu = nu/(enu-1) 
    :        = (z V F / (RT)) / (exp(z V F / (RT))-1)
    :        = (z V F / (RT))   (exp(-zV F / (RT))/(1-exp(-zV F / (RT))))
    :
    : now the effective drive is calculated as
    :
    :   ghkg = -f (1 - (ci/co)  enu) fnu
    :        = -10^3 RT/(zF)  (1 - (ci/co) exp(z V F / (RT))) *
    :         (z V F / (RT)) (exp(-zV F / (RT))/(1-exp(-zV F / (RT))))
    :        = -10^3 V (1/co) (co - ci exp(z V F / (RT))) (exp(-zV F / (RT))/(1-exp(-zV F / (RT))))
    :        = 10^3 V/co (ci - co exp(-zV F / (RT)))/(1-exp(-zV F / (RT)))
    :
    : [note, the 10^3 converts back to mV]
    : and you can see this is the ghk equation if the relationship
    : between conductance and permeability is
    :
    :      g = rho z^2 co F^2/RT
    :
    : Then g*ghkg reduces to the GHK current equation
    :
    
    ff   = (1.0e3/z)*R*(celsius+273.15)/FARADAY
    nu  = v/ff
	enu = exp(nu)
	if (fabs(nu) < 1e-4) {
        fnu = 1 - nu/2
    }else{
        fnu = nu/(enu-1) 
    }
    ghkg= -ff*(1-(ci/co)*enu)*fnu
}


