TITLE CA

:Includes changes made by Andreas Schaefer

:ca.mod 

:to lead to thalamic ca current inspired by destexhe and huguenrd                   
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX sit2
    USEION ca READ eca WRITE ica
    RANGE  gca, gamma, eta, deterministic,reff
    RANGE  N, inactF, actF  
    GLOBAL v12m, v12h, vwm, vwh, vm1, vm2, vh1, vh2, wm1, wm2, wh1, wh2 
    GLOBAL amc, ahc, minf, mtau, hinf, htau            
    GLOBAL vmin, vmax, vshift     
    GLOBAL P_am, P_bm, P_ah, P_bh, wflag
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
 
    v12m=50             (mV)
    v12h=78             (mV)
    vwm =7.4            (mV)
    vwh=5.0             (mV)
    amc=3            (mV)
    ahc=85           (mV)
    vm1=25          (mV)
    vm2=100             (mV)
    vh1=46          (mV)
    vh2=405             (mV)
    wm1=20          (mV)
    wm2=15          (mV)
    wh1=4           (mV)
    wh2=50          (mV)
                          
    gamma = 8      (pS)
    eta = 1      (1/um2)    : 0.0008 mho/cm2
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
    reff    (pS/um2)
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
    reff = eta*gamma
    N   = floor(eta*area + 0.5)   : round to nearest number of channels

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
    gca = m*m*h* N * scale_dens      
    } else { 
    gca = floor(m*m*h* N + 0.5) * scale_dens }
    } else{                                         
    gca = strap(m2h1) * scale_dens 
}
    ica = (1e-4) * gca * (v - eca) } 

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

PROCEDURE rates(v_(mV)) {     TABLE minf, hinf, mtau, htau, am, bm, ah, bh     
    DEPEND v12m,v12h,vwm,vwh,vm1,vm2,vh1,vh2,wm1,wm2,wh1,wh2,amc,ahc   
    FROM -120 TO 100 WITH 199

UNITSOFF
    minf = 1.0 / ( 1 + exp(-(v_+v12m)/vwm) )
    hinf = 1.0 / ( 1 + exp((v_+v12h)/vwh) )

    mtau = ( amc + 1.0 / ( exp((v_+vm1)/wm1) + exp(-(v_+vm2)/wm2) ) ) 
    htau = ( ahc + 1.0 / ( exp((v_+vh1)/wh1) + exp(-(v_+vh2)/wh2) ) )  
          
    ah = hinf/htau
    bh = 1/htau - ah                
    
    am = minf/mtau
    bm = 1/mtau - am   
      : no temperature adjustment done here
}

UNITSON  
: ----------------------------------------------------------------
: sign trap - trap negative numbers and replace with zero
FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"sit2.mod:strap: negative value for state");
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
    fprintf(stderr, "sit2.mod:ChkProb: argument not a probability.\n");
ENDVERBATIM
    wflag =0}
  }
}
