COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT


NEURON {

	SUFFIX cacurrent
	USEION ca WRITE ica
	RANGE gca, ica
}
	
UNITS {

	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gca = 10  (mS/cm2)
    eca = 80  (mV)
}
    
ASSIGNED { 

    ica  (mA/cm2)    
    v    (mV)
    sinf (1)
    taus (ms)
}

STATE { s }

INITIAL { 

    rates(v)
    s  = sinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ica = (1e-3) * gca * s^2 * (v-eca)
}


DERIVATIVE states { 

    rates(v)
    s' = (sinf-s)/taus
}

PROCEDURE rates(v(mV)) { LOCAL a,b

    a = fun2(v,      5,    1.6,     -1/0.072)
    b = fun3(v,   -8.9,   0.02,     5)

    sinf = a/(a+b)
    taus = 1.0/(a+b)
}

COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT



:-------------------------------------------------------------------
FUNCTION fun1(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun1 = A*exp((v-V0)/B)
}

FUNCTION fun2(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun2 = A/(exp((v-V0)/B)+1)
}

FUNCTION fun3(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

    if(fabs((v-V0)/B)<1e-6) {
    :if(v==V0) {
        fun3 = A*B/1(mV) * (1- 0.5 * (v-V0)/B)
    } else {
        fun3 = A/1(mV)*(v-V0)/(exp((v-V0)/B)-1)
    }
}

FUNCTION min(x,y) { if (x<=y){ min = x }else{ min = y } }