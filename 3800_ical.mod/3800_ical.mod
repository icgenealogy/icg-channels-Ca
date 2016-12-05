TITLE Cardiac L-type Calcium channel
: Hodgkin - Huxley type calcium channel from Courtemanche et al Am J Physiol 1998 275:H301b with voltage and calcium dependent inactivation

NEURON {
	SUFFIX ICaL
	USEION ca READ eca,cai WRITE ica
	RANGE gCaL, ica
	GLOBAL minf, ninf, hinf, mtau, ntau 
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (mM) = (milli/liter)
        
}

PARAMETER {
	gCaL=0.2476e-3 (S/cm2) <0,1e9> 
               
}

STATE {
	m n h
}

ASSIGNED {
        eca (mV)
	v (mV)
	celsius (degC) : 37
        cai (mM)
	ica (mA/cm2)
	minf ninf hinf
	mtau (ms)
	ntau (ms)
       
}

INITIAL {
	rates(v, cai)
	m = minf
	n = ninf
        h = hinf   
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	ica = gCaL*m*n*h*(v - eca)
}

DERIVATIVE states {
	rates(v, cai)
	m' = (minf - m)/mtau
        n' = (ninf - n)/ntau
	h' = (hinf - h)/2
}

FUNCTION alp(v(mV),i) (/ms) { LOCAL  q10
	q10 = 3^((celsius - 37(degC))/10(degC))
	if (i==0) {
		alp = (1 - exp(-(v + 10)/6.24))/(0.035*(v + 10)*(1 + exp(-(v + 10)/6.24)))/(q10*1(/ms))
	}else if (i==1){
		alp = 9/(0.0197*exp(-0.0337^2*(v + 10)^2) + 0.02)/(q10*1(/ms))
	}
}

FUNCTION bet(v(mV),i)(/ms) { 
	 v = v 
	
	if (i==0) {
		bet = 1/(1 + exp(-(v + 10)/8(mV)))
	}else if (i==1){
		bet = 1/(1 + exp((v + 28)/6.9(mV)))
	}
}

FUNCTION ce(cai(mM)) { 
            
           ce = 1/(1 + (cai/0.00035))
}



PROCEDURE rates(v(mV), cai (mM)) (/ms) {LOCAL a, b, c 
	
	a = alp(v,0)  b=bet(v,0)
	mtau = a
	minf = b
	a = alp(v,1)  b = bet(v,1)
	ntau = a
	ninf = b
        c = ce(cai)
	hinf = c
}
