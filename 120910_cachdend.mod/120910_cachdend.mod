TITLE cachdend.mod    
   
UNITS {  
        (mA) = (milliamp)  
        (mV) = (millivolt)  
}  
   
NEURON {  
        SUFFIX cachdend  
        USEION ca READ eca WRITE ica  
        RANGE gcabar 
        GLOBAL minf, mexp  
}  
   
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}  
   
PARAMETER {  
        v (mV)  
        celsius  (degC)  
        dt (ms)  
        gcabar = 4e-04 (mho/cm2)  
        eca = 125 (mV)  
}  
   
STATE {  
        m   
}  
   
ASSIGNED {  
        ica (mA/cm2)  
        minf mexp  
}  
   
BREAKPOINT {  
        SOLVE states  
        ica = gcabar*m*(v - eca)  
}  
   
UNITSOFF  
   
INITIAL {  
     rates(v)  
     m = minf  
}  
PROCEDURE states() {  :Computes state variable m  
        rates(v)      :             at the current v and dt.  
        m = m + mexp*(minf-m)  
}  
   
PROCEDURE rates(v) {:Computes rate and o  
         : ther constants at current v.  
         : Call once from HOC to   
         : initialize inf at resting v.  
     LOCAL  q10, tinc, alpha, beta, sum  
     TABLE minf, mexp DEPEND dt, celsius FROM -100 TO 100 WITH 200 

        q10 = 3^((celsius - 20)/10)  
        tinc = -dt * q10  
                :"m" calcium activation system  
        alpha = 1.5 * vtrap(-(v-20),5)  
        beta =  1.5 * exp(-(v+25)/10)  
        sum = alpha + beta  
        minf = alpha/sum  
        mexp = 1 - exp(tinc*sum)  
}  
   
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.  
        if (fabs(x/y) < 1e-6) {  
                vtrap = y*(1 - x/y/2)  
        }else{  
                vtrap = x/(exp(x/y) - 1)  
        }  
}  
   
UNITSON  
 