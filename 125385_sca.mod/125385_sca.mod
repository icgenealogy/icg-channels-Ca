TITLE 

:Includes changes made by Andreas Schaefer

:HVA Ca current
:Based on Reuveni, Friedman, Amitai and Gutnick (1993) J. Neurosci. 13:
:4609-4621.
                  
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX ssca
    USEION ca READ eca WRITE ica
    RANGE  gca, gamma, eta, deterministic 
    RANGE reff, N, inactF, actF
    GLOBAL  minf, mtau, hinf, htau     
    GLOBAL vmin, vmax, tadj 
    GLOBAL P_am, P_bm, P_ah, P_bh, wflag, vshift
    GLOBAL DONT_VECTORIZE   : prevent vectorization to agree with RNG.mod
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
    PI  = (pi) (1)   

}

PARAMETER {
    v (mV)
    dt (ms)
    area
    
    inactF = 3
    actF   = 1   

    gamma = 20      (pS)
    eta = .005      (1/um2)    : 
    vshift = 0  (mV)        : voltage shift (affects all)

    cao  = 2.5  (mM)            : external ca concentration
    cai     (mM)        celsius     (degC)
    vmin = -120 (mV)
    vmax = 100  (mV)

    DONT_VECTORIZE          : required declaration
    deterministic = 0   : if non-zero, use deterministic variables
}

STATE { m h                             : deterministic variables

m0h0 m0h1 m1h0 m1h1 m2h0 m2h1 : states

m0h0_m1h0  m1h0_m2h0  m0h1_m1h1  m1h1_m2h1  
m2h0_m1h0  m1h0_m0h0 m2h1_m1h1  m1h1_m0h1  
m0h0_m0h1 m0h1_m0h0 m1h0_m1h1 m1h1_m1h0 m2h0_m2h1 m2h1_m2h0
}

ASSIGNED {
    ica      (mA/cm2)
    gca      (pS/um2)
    eca      (mV)
    minf        
    hinf
    mtau (ms)   
    htau (ms)
    tadj

    am (/ms)
    bm (/ms)
    ah (/ms)
    bh (/ms)

    N 
    reff       (pS/um2)
    scale_dens (pS/um2)
    P_am        : probability of one channel making am transition
    P_bm
    P_ah
    P_bh

    wflag
}

INITIAL {
    rates(v+vshift)
    m = minf
    h = hinf 
    wflag = 1   : only give a warning once!

    scale_dens = gamma/area
    N   = floor(eta*area + 0.5)   : round to nearest number of channels
    reff = eta*gamma 

    m1h0 = floor(2*m*(1-m)*(1-h)*N + 0.5)
    m2h0 = floor(m*m*(1-h)*N + 0.5)

    m0h1 = floor((1-m)*(1-m)*h*N + 0.5)
    m1h1 = floor(2*m*(1-m)*h*N + 0.5)
    m2h1 = floor(m*m*h*N + 0.5)
    
    : put tbe rest of tbe channels in the non-conducting & inactivated state
    m0h0 = N - (m1h0 + m2h0 + m0h1 + m1h1 + m2h1)

    m0h0_m1h0=0 
    m1h0_m2h0=0 
    m0h1_m1h1=0
    m1h1_m2h1=0
    m2h0_m1h0=0
    m1h0_m0h0=0 
    m2h1_m1h1=0 
    m1h1_m0h1=0

    m0h0_m0h1=0 
    m0h1_m0h0=0 
    m1h0_m1h1=0 
    m1h1_m1h0=0 
    m2h0_m2h1=0 
    m2h1_m2h0=0 
}

BREAKPOINT {
  SOLVE states
    if (deterministic) { 
        if (deterministic-1){   
    gca = m*m*h*tadj*reff     
    } else { 
    gca = floor(m*m*h* N + 0.5) * scale_dens *tadj}
    } else{                                         
    gca = strap(m2h1) * scale_dens *tadj
}
    ica = (1e-4) * gca * (v - eca)
} 

PROCEDURE states() {
VERBATIM
    extern double BnlDev_RNG();
ENDVERBATIM
    rates(v+vshift)

    : deterministic versions of state variables
    : integrated by relaxing toward the steady state value
        m = m + (1 - exp(-dt/mtau)) * (minf-m)
        h = h + (1 - exp(-dt/htau)) * (hinf-h)    

    P_am = strap(am*dt)
    P_bm  = strap(bm*dt)
    
    : cbeck that will represent probabilities when used
    ChkProb( 2.0 * P_am)
    ChkProb( 2.0 * P_bm)
    ChkProb( P_bm/(1.0-P_am) )
    
    : m gate transitions

    m0h0_m1h0 = BnlDev_RNG(2.0*P_am,m0h0)
    m1h0_m2h0 = BnlDev_RNG(P_am,m1h0)
    m1h0_m0h0 = BnlDev_RNG(P_bm/(1.0-P_am), m1h0 - m1h0_m2h0)  
    m2h0_m1h0 = BnlDev_RNG(2.0*P_bm, m2h0)
    m0h1_m1h1 = BnlDev_RNG(2.0*P_am, m0h1)
    m1h1_m2h1 = BnlDev_RNG(P_am, m1h1)
    m1h1_m0h1 = BnlDev_RNG(P_bm/(1.0-P_am), m1h1 - m1h1_m2h1)
    m2h1_m1h1 = BnlDev_RNG(2.0*P_bm, m2h1)
    
    : new numbers in each state after the a gate transitions
    m0h0 = m0h0 - m0h0_m1h0 + m1h0_m0h0
    m1h0 = m1h0 - m1h0_m2h0 - m1h0_m0h0  + m2h0_m1h0 + m0h0_m1h0
    m2h0 = m2h0 - m2h0_m1h0 + m1h0_m2h0
 
    m0h1 = m0h1 - m0h1_m1h1 + m1h1_m0h1
    m1h1 = m1h1 - m1h1_m2h1 - m1h1_m0h1 + m2h1_m1h1 + m0h1_m1h1
    m2h1 = m2h1 - m2h1_m1h1 + m1h1_m2h1

    : probabilities of making h gate transitions
    P_ah = strap(ah*dt)
    P_bh  = strap(bh*dt)
    
    ChkProb(P_ah)
    ChkProb(P_bh)
    
    : number making h gate transitions

    m0h0_m0h1 = BnlDev_RNG(P_ah,m0h0)
    m0h1_m0h0 = BnlDev_RNG(P_bh,m0h1)
    m1h0_m1h1 = BnlDev_RNG(P_ah,m1h0)
    m1h1_m1h0 = BnlDev_RNG(P_bh,m1h1)
    m2h0_m2h1 = BnlDev_RNG(P_ah,m2h0)
    m2h1_m2h0 = BnlDev_RNG(P_bh,m2h1)

    m0h0 = m0h0 - m0h0_m0h1  + m0h1_m0h0
    m1h0 = m1h0 - m1h0_m1h1  + m1h1_m1h0
    m2h0 = m2h0 - m2h0_m2h1  + m2h1_m2h0

    m0h1 = m0h1 - m0h1_m0h0  + m0h0_m0h1
    m1h1 = m1h1 - m1h1_m1h0  + m1h0_m1h1
    m2h1 = m2h1 - m2h1_m2h0  + m2h0_m2h1

}

PROCEDURE rates(vm(mV)) {     
    TABLE minf, hinf, mtau, htau, am, bm, ah, bh, tadj
    DEPEND celsius, actF, inactF  
    FROM -120 TO 100 WITH 199

UNITSOFF
        tadj = 2.3^((celsius - 23)/10)

    am = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1)/actF
    bm = 0.94*exp((-75-vm)/17)/actF
    am = am*tadj
    bm = bm*tadj
    
    mtau = 1/(am+bm)     
    minf = am/(am+bm)

        :"h" inactivation 

    ah = 0.000457*exp((-13-vm)/50)/inactF
    bh = 0.0065/(exp((-vm-15)/28) + 1)/inactF
    ah = ah*tadj
    bh = bh*tadj

    htau = 1/(ah+bh)
    hinf = ah/(ah+bh)

}

UNITSON  
: ----------------------------------------------------------------
: sign trap - trap negative numbers and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"sca.mod:strap: negative value for state");
ENDVERBATIM
    } else { 
        strap = x
    }
}

: ----------------------------------------------------------------
: ChkProb - Check that number represents a probability
PROCEDURE ChkProb(p) {
  if (p < 0.0 || p > 1.0) {
    if (wflag){
VERBATIM
    fprintf(stderr, "sca.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
    wflag =0}
  }
}
