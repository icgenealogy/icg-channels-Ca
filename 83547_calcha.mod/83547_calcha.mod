TITLE calcha.mod 
 
COMMENT
ENDCOMMENT
 
UNITS {
       (molar) = (1/liter)
        (S) = (siemens)
         (mA) = (milliamp)
        (mV) = (millivolt)
         F = (faraday) (coulomb)
        R = (mole k) (mV-coulomb/degC)
       (mM) =  (millimolar)

}
 
NEURON {
        SUFFIX calcha
        USEION ca READ eca, cai WRITE ica
        RANGE gcalbar, gcanbar, gcahvabar, gcatbar,ica, ical,icat,ican,icahva,kml,kmn
        GLOBAL dlinf, dninf,  dtinf,ftinf,dhvainf,fhvainf
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 35.0 (degC)
        dt (ms)
        gcatbar =  1044.0e-6 (S/cm2)
        gcanbar =  171.0e-6  (S/cm2)
        gcalbar =  216.0e-6  (S/cm2)
        gcahvabar =  0.0e-6  (S/cm2)
        kmn = 0.0001   (mM)
        kml = 0.00045  (mM)
        eca (mV)
        cao = 2.0 (mM)
        cai   (mM)

}
 
STATE {
        d_t dl dn dhva fhva ft
}
 
ASSIGNED {
        ica (mA/cm2)
        icahva (mA/cm2)
        ical (mA/cm2)
        ican (mA/cm2)
        icat (mA/cm2)
        dlinf dninf dtinf ftinf dhvainf fhvainf fninf flinf
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        fninf = kmn/(kmn + cai)
        flinf = kml/(kml + cai)
        COMMENT
        :eca = R*(celsius+273.15)/F*log(cao/cai)
        ENDCOMMENT
        icahva = gcahvabar*dhva*fhva*(v - eca)
        ical = gcalbar*dl*flinf*(v - eca)
        ican = gcanbar*dn*fninf*(v - eca)
        icat = gcatbar*d_t*ft*(v - eca)
        ica = icahva + ical + ican + icat
}
 
UNITSOFF
 
INITIAL {
        dhva = boltz(v,-10.0,10.0)
        dl = boltz(v,-50.0,3.0)
        dn = boltz(v,-45.0,7.0)
        d_t = boltz(v,-63.5,1.5)
        fhva = boltz(v,-48.0,-5.0)
        ft = boltz(v,-76.2,-3.0)
}

DERIVATIVE states {  :Computes state variables m, h, and n 
LOCAL dlinf,dninf,dtinf,ftinf,dhvainf,fhvain,dltau,dntau,dttau,fttau,dhvatau,fhvatau
        dhvainf = boltz(v,-10.0,10.0)
        dlinf = boltz(v,-50.0,3.0)
        dninf = boltz(v,-45.0,7.0)
        dtinf = boltz(v,-63.5,1.5)
        fhvainf = boltz(v,-48.0,-5.0)
        ftinf = boltz(v,-76.2,-3.0)
        dhvatau = gaussian(v,0.1,13.0,62.0,0.05)
        dltau = gaussian(v,18.0,20.0,45.0,1.50)
        dntau = gaussian(v,18.0,25.0,70.0,0.30)
        dttau = gaussian(v,65.0,6.32455,66.0,3.5)
        fhvatau = gaussian(v,0.5,18.0,55.6,0.50)
        fttau = gaussian(v,50.0,10.0,72.0,10.0)
        dhva' = (dhvainf-dhva)/dhvatau
        dl' = (dlinf-dl)/dltau
        dn' = (dninf-dn)/dntau
        d_t' = (dtinf-d_t)/dttau
        fhva' = (fhvainf-fhva)/fhvatau
        ft' = (ftinf-ft)/fttau
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

