TITLE calcium channels (L, N, and T types) 
 
UNITS {
       (molar) = (1/liter)
       (S)  = (siemens)
       (mA) = (milliamp)
       (mV) = (millivolt)
       (mM) = (millimolar)
        F = (faraday)  (coulomb)
        R = (mole k)   (mV-coulomb/degC)
       
}
 
NEURON {
        SUFFIX cachan
        USEION ca READ eca,cai WRITE ica
        RANGE  gcalbar,gcanbar,gcatbar,ica,ical,icat,ican,kml,kmn
        GLOBAL dlinf,dninf,dtinf,ftinf,flinf,fninf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        cai   (mM)
        celsius =  35.0      (degC)
        gcatbar =  1044.0e-6 (S/cm2)
        gcanbar =  171.0e-6  (S/cm2)
        gcalbar =  216.0e-6  (S/cm2)
        kmn = 0.0001   (mM)
        kml = 0.00045  (mM)
        eca (mV)
        cao = 2.0 (mM)
        
}
 
STATE {
        d_t dl dn ft
}
 
ASSIGNED {
        ica (mA/cm2)
        ical (mA/cm2)
        ican (mA/cm2)
        icat (mA/cm2)
        dlinf dninf dtinf ftinf fninf flinf
 }
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        fninf = kmn/(kmn + cai)
        flinf = kml/(kml + cai)
        ical = gcalbar*dl*flinf*(v - eca)
        ican = gcanbar*dn*fninf*(v - eca)
        icat = gcatbar*d_t*ft*(v - eca)
        ica  = ical + ican + icat
}
 
UNITSOFF
 
INITIAL {
        dl = boltz(v,-50.0,3.0)
        dn = boltz(v,-45.0,7.0)
        d_t = boltz(v,-63.5,1.5)
        ft = boltz(v,-76.2,-3.0)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL dlinf,dninf,dtinf,ftinf,dltau,dntau,dttau,fttau
        dlinf = boltz(v,-50.0,3.0)
        dninf = boltz(v,-45.0,7.0)
        dtinf = boltz(v,-63.5,1.5)
        ftinf = boltz(v,-76.2,-3.0)
        dltau = gaussian(v,18.0,20.0,45.0,1.50)
        dntau = gaussian(v,18.0,25.0,70.0,0.30)
        dttau = gaussian(v,65.0,6.32455,66.0,3.5)
        fttau = gaussian(v,50.0,10.0,72.0,10.0)
        dl'  = (dlinf-dl)/dltau
        dn'  = (dninf-dn)/dntau
        d_t' = (dtinf-d_t)/dttau
        ft'  = (ftinf-ft)/fttau
}
 
 
FUNCTION gaussian(v,a,b,c,d) {
        LOCAL arg
        arg= a*exp(-(c+v)*(v+c)/(b*b)) +d
        gaussian = arg
}
 
 
FUNCTION boltz(x,y,z) {
               LOCAL arg
                arg= -(x-y)/z
                if (arg > 50) {boltz = 0}
                else {if (arg < -50) {boltz = 1}
                else {boltz = 1.0/(1.0 + exp(arg))}}
}

 
UNITSON

