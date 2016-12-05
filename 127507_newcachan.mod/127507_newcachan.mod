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
        USEION ca READ eca WRITE ica
        RANGE  gcalbar,gcanbar,gcahvabar,ica,ical,icahva,ican,kml,kmn
        GLOBAL dlinf,dhvainf,fhvainf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        dt (ms)
        :cai   (mM)
        celsius =  35.0      (degC)
        gcahvabar =  0.0e-6 (S/cm2)
        gcalbar =  11.196e-6  (S/cm2)
        kmn = 0.0001   (mM)
        kml = 0.00045  (mM)
        :eca = 120 (mV)
        cao = 2.0 (mM)
        
}
 
STATE {
        dhva dl  fhva
}
 
ASSIGNED {
        eca (mV)
        ica (mA/cm2)
        ical (mA/cm2)
        ican (mA/cm2)
        icahva (mA/cm2)
        dlinf  dhvainf fhvainf 
 }
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        ical = gcalbar*dl*(v - eca)
        icahva = gcahvabar*dhva*fhva*(v - eca)
        ica  = ical +  icahva
}
 
UNITSOFF
 
INITIAL {
        dl = boltz(v,-45.0,5.0)
        dhva = boltz(v,-10.0,10.0)
        fhva = boltz(v,-48.0,-5.0)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL dlinf,dhvainf,fhvainf,dltau,dhvatau,fhvatau
        dlinf = boltz(v,-45.0,5.0)
        dhvainf = boltz(v,-10.0,10.0)
        fhvainf = boltz(v,-48.0,-5.0)
        dltau = gaussian(v,18.0,25.0,70.0,0.30)
        dhvatau = gaussian(v,0.1,13.0,62.0,0.05)
        fhvatau = gaussian(v,0.5,18.0,55.6,0.5)
        dl'  = (dlinf-dl)/dltau
        dhva' = (dhvainf-dhva)/dhvatau
        fhva'  = (fhvainf-fhva)/fhvatau
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

